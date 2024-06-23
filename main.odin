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
    return unit(float3{rand.float32_range(-1.0, 1.0), rand.float32_range(-1.0, 1.0), rand.float32_range(-1.0, 1.0)})
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

near_zero :: proc(v: float3) -> bool {
    ep : f32 = 1e-8
    return (math.abs(v.x) <= ep) || (math.abs(v.y) <= ep) || (math.abs(v.z) <= ep)
}

reflect :: proc(v: float3, n: float3) -> float3 {
    return v - float3(2.0) * dot(v, n) * n
}

//Attempts to return the refracted ray, returns the reflected ray
//if refraction isn't possible
refract :: proc(v: float3, n: float3, ri: f32) -> float3 {
    cos_theta := math.min(dot(-v, n), 1.0)
    sin_theta := math.sqrt(1.0 - cos_theta*cos_theta)

    if (ri * sin_theta > 1.0 || schlick_reflectance(cos_theta, ri) > rand.float32()) {
        //Must reflect
        return reflect(v, n)
    } else {
        //Can refract
        refract_perpendicular := float3(ri) * (v + cos_theta * n)
        refract_parallel := float3(-math.sqrt(1.0 - length_squared(refract_perpendicular))) * n
        return refract_perpendicular + refract_parallel
    }
}

schlick_reflectance :: proc(cosine: f32, ri: f32) -> f32 {
    one_min_cos := (1.0 - cosine)
    r0 := (1.0 - ri) / (1.0 + ri)
    r0 = r0*r0
    return r0 + (1.0 - r0) * one_min_cos * one_min_cos * one_min_cos * one_min_cos * one_min_cos
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
    hit_sphere_idx : u64,
    material_idx: u64,
    front_face: bool
}

sphere :: struct {
    origin : float3,
    radius : f32,
    material_idx: u64
}

material_type :: enum {Dielectric, Matte, Metal }
material :: struct {
    albedo: float3,
    roughness: f32,
    refractive_index: f32,
    type: material_type
}

scene :: struct {
    spheres: [dynamic]sphere,
    materials: [dynamic]material
}

linear_to_gamma :: proc(color: float3) -> float3 {
    return float3 {
        math.sqrt(color.r),
        math.sqrt(color.g),
        math.sqrt(color.b),
    }
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
    origin_difference := s.origin - r.origin
    a := length_squared(r.direction)
    h := dot(r.direction, origin_difference)
    c := length_squared(origin_difference) - s.radius * s.radius
    discriminant := h*h - a*c

    if discriminant < 0 do return nil

    //We only do - instead of + as well because
    //we are looking for the closest intersection point
    t := (h - math.sqrt(discriminant)) / a
    if t < min_t {
        //Unless we could be inside the sphere
        t = (h + math.sqrt(discriminant)) / a
        if t < min_t do return nil
    }
    return t
}

sample_square :: proc() -> float3 {
    x := rand.float32()
    y := rand.float32()
    return float3{x - 0.5, y - 0.5, 0.0}
}

ray_color :: proc(r: ray, bounces_left: u8, scene: scene) -> float3 {
    if bounces_left == 0 do return float3(0.0)

    //Trace against all spheres, updating the payload as we go
    payload := ray_payload {
        t_val = math.INF_F32
    }
    for sphere, idx in scene.spheres {
        res := ray_hit_sphere(r, sphere, 0.001)
        if res != nil && res.? < payload.t_val {
            payload.t_val = res.?
            payload.hit_sphere_idx = u64(idx)
            payload.material_idx = sphere.material_idx
        }
    }

    if payload.t_val < math.INF_F32 {
        //Now payload has closest-intersection info
        //Can think of this as the shading step
        intersection_point := r.origin + float3(payload.t_val) * r.direction
        outward_normal := unit(intersection_point - scene.spheres[payload.hit_sphere_idx].origin)
        payload.front_face = dot(r.direction, outward_normal) < 0
        
        material := scene.materials[payload.material_idx]
        new_ray: ray
        attenuation: float3
        
        switch material.type {
            case .Matte:
                attenuation = material.albedo
                new_ray_direction := outward_normal + rand_float3_sphere()
                if near_zero(new_ray_direction) do new_ray_direction = outward_normal
        
                new_ray = ray {
                    origin = intersection_point,
                    direction = new_ray_direction
                }
            case .Metal:
                attenuation = material.albedo
                new_ray_direction := reflect(r.direction, outward_normal)
                new_ray_direction += float3(material.roughness) * rand_float3_sphere()

                new_ray = ray {
                    origin = intersection_point,
                    direction = new_ray_direction
                }
            case .Dielectric:
                attenuation = float3(1.0)
                ri := payload.front_face ? (1.0 / material.refractive_index) : material.refractive_index
                outward_normal = payload.front_face ? outward_normal : -outward_normal
                refracted_ray := refract(unit(r.direction), outward_normal, ri)

                new_ray = ray {
                    origin = intersection_point,
                    direction = refracted_ray
                }
        }
        return attenuation * ray_color(new_ray, bounces_left - 1, scene)
    } else {
        //Ray missed all objects so return sky color
        unit_ray_dir := unit(r.direction)
        a := float3(0.5 * (unit_ray_dir.y + 1.0))
        return (float3(1.0) - a) * float3{1.0, 1.0, 1.0} + a * float3{0.5, 0.7, 1.0}
    }

    return float3(0.0)
}

main :: proc() {
    fmt.println("Initiating swag mode...")

    image_width := 640
    image_height := 360
    aspect_ratio := f32(image_width) / f32(image_height)
    samples_per_pixel := 100
    max_ray_depth : u8 = 10

    //Internal framebuffer we hold in RAM during rendering
    framebuffer : [dynamic]float3
    resize(&framebuffer, image_height * image_width)

    height : f32 = 2.0
    camera := camera {
        origin = float3{0.0, 0.0, 0.0},
        focal_length = 1.0,
        viewport_height = height,
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
    scene := scene {
        spheres = [dynamic]sphere{
            {
                origin = float3{0.0, 0.0, -1.2},
                radius = 0.5,
                material_idx = 0
            },
            {
                origin = float3{1.0, 0.0, -1.0},
                radius = 0.5,
                material_idx = 2
            },
            {
                origin = float3{-1.0, 0.0, -1.0},
                radius = 0.5,
                material_idx = 3
            },
            {
                origin = float3{-1.0, 0.0, -1.0},
                radius = 0.4,
                material_idx = 5
            },
            {
                origin = float3{0, -100.5, -1.0},
                radius = 100.0,
                material_idx = 4
            }
        },
        materials = [dynamic]material {
            {
                albedo = float3{0.1, 0.2, 0.5},
                roughness = 1.0,
                type = .Matte
            },
            {
                albedo = float3(0.8),
                roughness = 0.4,
                type = .Metal
            },
            {
                albedo = float3{0.8, 0.6, 0.2},
                roughness = 1.0,
                type = .Metal
            },
            {
                refractive_index = 1.5,
                type = .Dielectric
            },
            {
                albedo = float3{0.8, 0.8, 0.0},
                type = .Matte
            },
            {
                refractive_index = 1.0 / 1.5,
                type = .Dielectric
            }
        }
    }

    //Iterate through the image's pixels
    //from top-left to bottom-right, 
    fmt.println("Beginning main tracing routine...")
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
                ray_direction := pixel_center - camera.origin
                ray := ray {
                    origin = camera.origin,
                    direction = ray_direction
                }
    
                color += ray_color(ray, max_ray_depth, scene)
            }
            color /= float3(samples_per_pixel)

            //Write color to correct pixel
            pixel_idx := j * image_width + i
            framebuffer[pixel_idx] = color
        }
    }

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
            pixel = linear_to_gamma(pixel)
            red_int := i32(pixel.r * 255.0)
            green_int := i32(pixel.g * 255.0)
            blue_int := i32(pixel.b * 255.0)
            fmt.fprintln(image_file, red_int, green_int, blue_int, flush = false)
        }
    }

    fmt.println("Done!")
}