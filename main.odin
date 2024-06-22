//Author: Nick Driscoll
//Descripton: Implementation of Ray Tracing in One Weekend
//https://raytracing.github.io/books/RayTracingInOneWeekend.html

package main

import "core:math"
import "core:fmt"
import "core:os"
import "vendor:stb/image"

float3 :: distinct [3]f32
dot :: proc(a: float3, b: float3) -> f32 {
    return a.x * b.x + a.y * b.y + a.z * b.z
}

ray :: struct {
    origin : float3,
    direction : float3
}
unit :: proc(vec : float3) -> float3 {
    term := vec.x * vec.x + vec.y * vec.y + vec.z * vec.z
    length := math.sqrt(term)
    return vec / float3(length)
}

ray_payload :: struct {
    t_val : f32,
    color : float3,
    hit_sphere_idx : u64
}

sphere :: struct {
    //color : float3
    origin : float3,
    radius : f32
}

//Returns the t-value for closest intersection
ray_hit_sphere :: proc(r: ray, s: sphere) -> Maybe(f32) {
    //Directly solve the quadratic
    origin_difference := s.origin - r.origin
    a := dot(r.direction, r.direction)
    b := -2.0 * dot(r.direction, origin_difference)
    c := dot(origin_difference, origin_difference) - s.radius * s.radius

    sqrt_term := b * b - 4 * a * c
    if sqrt_term < 0 do return nil      //Early exit if the quadratic has no real solutions

    //We only do - instead of + as well because
    //we are looking for the closest intersection point
    t := (-b - math.sqrt(sqrt_term)) / (2 * a)
    if t < 0 do return nil
    return t
}

main :: proc() {
    fmt.println("Initiating swag mode...")
    ctxt := context

    image_width := 640
    image_height := 480
    aspect_ratio := f32(image_width) / f32(image_height)

    //Internal framebuffer we hold in RAM during rendering
    framebuffer : [dynamic]float3
    resize(&framebuffer, image_height * image_width)

    focal_length : f32 = 1.0
    viewport_height : f32 = 2.0
    viewport_width := aspect_ratio * viewport_height
    camera_origin := float3{0.0, 0.0, 0.0}

    //Horizontal and vertical viewport unit vectors
    viewport_u := float3{viewport_width, 0.0, 0.0}
    viewport_v := float3{0.0, -viewport_height, 0.0}
    pixel_delta_u := viewport_u / float3(image_width)
    pixel_delta_v := viewport_v / float3(image_height)
    
    viewport_top_left := camera_origin - float3{0.0, 0.0, focal_length} - viewport_u / 2.0 - viewport_v / 2.0
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
        fmt.print((image_height - j),"lines remaining...     \r")
        for i := 0; i < image_width; i += 1 {
            //Compute color
            color : float3

            pixel_center := pixel00_center + (float3(i) * pixel_delta_u) + (float3(j) * pixel_delta_v)
            ray_direction := pixel_center - camera_origin;
            ray := ray {
                origin = camera_origin,
                direction = ray_direction
            }

            payload := ray_payload {
                t_val = math.INF_F32,
                color = float3(0.0)
            }
            for sphere, idx in spheres {
                res := ray_hit_sphere(ray, sphere)
                if res != nil && res.? < payload.t_val {
                    payload.t_val = res.?
                    payload.hit_sphere_idx = u64(idx)
                }
            }

            //Now payload has closest-intersection info
            if payload.t_val < math.INF_F32 {
                intersection_point := ray.origin + float3(payload.t_val) * ray.direction
                sphere_normal := unit(intersection_point - spheres[payload.hit_sphere_idx].origin)
                color = float3(0.5) * sphere_normal + float3(0.5)
            } else {
                //Sky color
                unit_ray_dir := unit(ray_direction)
                a := float3(0.5 * (unit_ray_dir.y + 1.0))
                color = (float3(1.0) - a) * float3{1.0, 1.0, 1.0} + a * float3{0.5, 0.7, 1.0}
            }

            // res := ray_hit_sphere(ray, main_sphere)
            // if res != nil {
            //     intersection_point := ray.origin + float3(res.?) * ray.direction
            //     sphere_normal := unit(intersection_point - main_sphere.origin)
            //     color = float3(0.5) * sphere_normal + float3(0.5)
            // } else {
            //     //Sky color
            //     unit_ray_dir := unit(ray_direction)
            //     a := float3(0.5 * (unit_ray_dir.y + 1.0))
            //     color = (float3(1.0) - a) * float3{1.0, 1.0, 1.0} + a * float3{0.5, 0.7, 1.0}
            // }

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