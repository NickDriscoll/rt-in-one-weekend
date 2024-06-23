//Author: Nick Driscoll
//Descripton: Implementation of Ray Tracing in One Weekend
//https://raytracing.github.io/books/RayTracingInOneWeekend.html

package main

import "core:math"
import "core:math/rand"
import "core:fmt"
import "core:os"
import "vendor:stb/image"

float3 :: distinct [3]f32
dot :: proc(a: float3, b: float3) -> f32 {
    return a.x * b.x + a.y * b.y + a.z * b.z
}
length_squared :: proc(vec: float3) -> f32 {
    return dot(vec, vec)
}
length :: proc(vec: float3) -> f32 {
    return math.sqrt(length_squared(vec))
}
unit :: proc(vec : float3) -> float3 {
    l := length(vec)
    return vec / float3(l)
}

rand_float3_sphere :: proc() -> float3 {
    return float3{rand.float32_range(-1.0, 1.0), rand.float32_range(-1.0, 1.0), rand.float32_range(-1.0, 1.0)}
}
rand_unit_float3 :: proc() -> float3 {
    return unit(rand_float3_sphere())
}

//Returns a random unit vector in the hemisphere defined by normal
rand_float3_hemisphere :: proc(normal: float3) -> float3 {
    v := rand_float3_sphere()
    v = unit(v)
    if dot(v, normal) < 0 {
        v = -v
    }
    return v
}

interval :: struct {
    min: f32,
    max: f32
}

camera :: struct {
    origin: float3,
    focal_length: f32,
    viewport_width: f32,
    viewport_height: f32
}

ray :: struct {
    origin : float3,
    direction : float3
}

ray_payload :: struct {
    t_val : f32,
    hit_sphere_idx : u64
}

sphere :: struct {
    origin : float3,
    radius : f32
}

//Returns the t-value for closest intersection
ray_hit_sphere_naive :: proc(r: ray, s: sphere) -> Maybe(f32) {
    //Directly solve the quadratic
    origin_difference := s.origin - r.origin
    a := dot(r.direction, r.direction)
    b := -2.0 * dot(r.direction, origin_difference)
    c := dot(origin_difference, origin_difference) - s.radius * s.radius

    sqrt_term := b * b - 4 * a * c
    if sqrt_term < 0.0 do return nil      //Early exit if the quadratic has no real solutions

    //We only do - instead of + as well because
    //we are looking for the closest intersection point
    t := (-b - math.sqrt(sqrt_term)) / (2 * a)
    if t < 0.0 do return nil
    return t
}

ray_hit_sphere :: proc(r: ray, s: sphere, min_t : f32 = 0.0) -> Maybe(f32) {
    origin_difference := s.origin - r.origin;
    a := length_squared(r.direction)
    h := dot(r.direction, origin_difference)
    c := length_squared(origin_difference) - s.radius * s.radius
    discriminant := h*h - a*c

    if discriminant < 0 do return nil

    //We only do - instead of + as well because
    //we are looking for the closest intersection point
    t := (h - math.sqrt(discriminant)) / a
    if t < min_t do return nil
    return t
}

sample_square :: proc() -> float3 {
    x := rand.float32()
    y := rand.float32()
    return float3{x - 0.5, y - 0.5, 0.0}
}

ray_color :: proc(r: ray, depth: u8, spheres: []sphere) -> float3 {
    if depth == 0 do return float3(0.0)

    //Trace against all spheres, updating the payload as we go
    payload := ray_payload {
        t_val = math.INF_F32
    }
    for sphere, idx in spheres {
        res := ray_hit_sphere(r, sphere, 0.001)
        if res != nil && res.? < payload.t_val {
            payload.t_val = res.?
            payload.hit_sphere_idx = u64(idx)
        }
    }

    //Now payload has closest-intersection info
    //Can think of this as the shading step
    if payload.t_val < math.INF_F32 {
        intersection_point := r.origin + float3(payload.t_val) * r.direction
        sphere_normal := unit(intersection_point - spheres[payload.hit_sphere_idx].origin)
        
        new_ray := ray {
            origin = intersection_point,
            direction = sphere_normal + rand_unit_float3()
        }
        return 0.5 * ray_color(new_ray, depth - 1, spheres)


        //color = float3(0.5) * sphere_normal + float3(0.5) //Just returning the normal vector as a color
    } else {
        //Ray missed all objects so return sky color
        unit_ray_dir := unit(r.direction)
        a := float3(0.5 * (unit_ray_dir.y + 1.0))
        return (float3(1.0) - a) * float3{1.0, 1.0, 1.0} + a * float3{0.5, 0.7, 1.0}
    }
}

main :: proc() {
    fmt.println("Initiating swag mode...")

    image_width := 640
    image_height := 480
    aspect_ratio := f32(image_width) / f32(image_height)
    samples_per_pixel := 10
    max_ray_depth : u8 = 10

    //Internal framebuffer we hold in RAM during rendering
    framebuffer : [dynamic]float3
    resize(&framebuffer, image_height * image_width)

    height : f32 = 2.0
    camera := camera {
        origin = float3{0.0, 0.0, 0.0},
        focal_length = 1.0,
        viewport_height = 2.0,
        viewport_width = aspect_ratio * height
    }

    //Horizontal and vertical viewport unit vectors
    viewport_u := float3{camera.viewport_width, 0.0, 0.0}
    viewport_v := float3{0.0, -camera.viewport_height, 0.0}
    pixel_delta_u := viewport_u / float3(image_width)
    pixel_delta_v := viewport_v / float3(image_height)
    
    viewport_top_left := camera.origin - float3{0.0, 0.0, camera.focal_length} - viewport_u / 2.0 - viewport_v / 2.0
    pixel00_center := viewport_top_left + 0.5 * (pixel_delta_u + pixel_delta_v)

    //Define scene objects
    spheres := [?]sphere{
        {
            origin = float3{0.0, 0.0, -1.0},
            radius = 0.5
        },
        {
            origin = float3{0, -100.5, -1.0},
            radius = 100.0
        }
    }
    main_sphere := sphere {
        origin = float3{0.0, 0.0, -1.0},
        radius = 0.5
    }

    //Iterate through the image's pixels
    //from top-left to bottom-right, 
    for j := 0; j < image_height; j += 1 {
        fmt.print((image_height - j), "lines remaining...     \r")
        for i := 0; i < image_width; i += 1 {
            //Compute color
            color := float3(0.0)

            //We send multiple, randomly jittered rays through the image plane
            for s := 0; s < samples_per_pixel; s += 1 {
                //Generate the ray to send through the current pixel
                jitter_vec := sample_square()
                pixel_center := pixel00_center + ((float3(i) + jitter_vec.x) * pixel_delta_u) + ((float3(j) + jitter_vec.y) * pixel_delta_v)
                ray_direction := pixel_center - camera.origin;
                ray := ray {
                    origin = camera.origin,
                    direction = ray_direction
                }
    
                color += ray_color(ray, max_ray_depth, spheres[:])
            }
            color /= float3(samples_per_pixel)

            //Write color to correct pixel
            pixel_idx := j * image_width + i
            framebuffer[pixel_idx] = color
        }
    }
    fmt.println()

    //File writing time, baby
    image_file, err := os.open("raytraced_image.ppm", os.O_CREATE | os.O_TRUNC)
    if err != os.ERROR_NONE {
        fmt.println("There was an error opening the file:", err)
    }
    defer os.close(image_file)

    //Write PPM header
    os.write_string(image_file, "P3\n")
    fmt.fprintln(image_file, image_width, image_height)
    os.write_string(image_file, "255\n")

    fmt.println("Writing framebuffer to file...")
    for j := 0; j < image_height; j += 1 {
        for i := 0; i < image_width; i += 1 {
            //One can imagine SIMDing the hell out of this
            pixel_idx := j * image_width + i
            pixel := framebuffer[pixel_idx]
            red_int := i32(pixel.r * 255.0)
            green_int := i32(pixel.g * 255.0)
            blue_int := i32(pixel.b * 255.0)
            fmt.fprintln(image_file, red_int, green_int, blue_int, flush = false)
        }
    }

    fmt.println("Done!")
}