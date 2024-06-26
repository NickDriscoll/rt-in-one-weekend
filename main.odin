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
cross :: proc(a: float3, b: float3) -> float3 {
    return float3 {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    }
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

rand_float3 :: proc() -> float3 {
    return float3{rand.float32_range(-1.0, 1.0), rand.float32_range(-1.0, 1.0), rand.float32_range(-1.0, 1.0)}
}

rand_float3_range :: proc(min: f32, max: f32) -> float3 {
    return float3{rand.float32_range(min, max), rand.float32_range(min, max), rand.float32_range(min, max)}
}

rand_float3_disk :: proc() -> float3 {
    for {
        p := float3{rand.float32_range(-1.0, 1.0), rand.float32_range(-1.0, 1.0), 0.0}
        if (length_squared(p) <= 1.0) {
            return p;
        }
    }
}

rand_float3_sphere :: proc() -> float3 {
    return unit(rand_float3())
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
    pixel_delta_u: float3,
    pixel_delta_v: float3,
    pixel00_center: float3,
    defocus_disk_u: float3,
    defocus_disk_v: float3,
    u, v, w : float3,       //Camera frame orthonormal basis vectors
    lookat: float3,         //Point camera is looking towards
    focal_length: f32,
    image_width: int,
    image_height: int,
    viewport_width: f32,
    viewport_height: f32,
    vertical_fov: f32,
    framebuffer: [dynamic]float3
}
init_camera :: proc(image_x: int, image_y: int, origin: float3, lookat: float3, v_fov: f32, defocus_angle: f32) -> camera {
    h : f32 = math.tan(v_fov / 2.0)

    focal_length := length(origin - lookat)
    viewport_height := 2.0 * h * focal_length
    viewport_width := (f32(image_x) / f32(image_y)) * viewport_height
    
    //Calculate basis vectors
    vup := float3{0.0, 1.0, 0.0}
    w := unit(origin - lookat)
    u := unit(cross(vup, w))
    v := cross(w, u)
    
    //Horizontal and vertical viewport unit vectors
    viewport_u := viewport_width * u
    viewport_v := -viewport_height * v
    pixel_delta_u := viewport_u / float3(image_x)
    pixel_delta_v := viewport_v / float3(image_y)

    viewport_top_left := origin - (focal_length * w) - viewport_u / 2.0 - viewport_v / 2.0
    pixel00_center := viewport_top_left + 0.5 * (pixel_delta_u + pixel_delta_v)

    defocus_radius := focal_length * math.tan(defocus_angle / 2.0);
    defocus_disk_u := u * defocus_radius
    defocus_disk_v := v * defocus_radius

    //Internal framebuffer we hold in RAM during rendering
    framebuffer : [dynamic]float3
    resize(&framebuffer, camera.image_height * camera.image_width)

    return camera {
        origin = origin,
        pixel_delta_u = pixel_delta_u,
        pixel_delta_v = pixel_delta_v,
        pixel00_center = pixel00_center,
        defocus_disk_u = defocus_disk_u,
        defocus_disk_v = defocus_disk_v,
        u = u, v = v, w = w,
        lookat = lookat,
        focal_length = focal_length,
        image_width = image_x,
        image_height = image_y,
        viewport_height = viewport_height,
        viewport_width = viewport_width,
        vertical_fov = v_fov,
        framebuffer = framebuffer
    }
}
get_pixel_ray :: proc(cam: camera, x_coord: u32, y_coord: u32) -> ray {
    //Generate the ray to send through the current pixel
    jitter_vec := sample_square()
    pixel_center := cam.pixel00_center + ((float3(x_coord) + jitter_vec.x) * cam.pixel_delta_u) + ((float3(y_coord) + jitter_vec.y) * cam.pixel_delta_v)
    random_unit_disk := rand_float3_disk()
    ray_origin := cam.origin + cam.defocus_disk_u * random_unit_disk.x + cam.defocus_disk_v * random_unit_disk.y
    ray_direction := pixel_center - ray_origin
    return ray {
        origin = ray_origin,
        direction = ray_direction
    }
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

//Returns the t-value for the closest intersection point on the sphere, or nil if it missed
//I.E the value 't' in the expression "intersection_point = r.origin + t * r.direction"
//This procedure accounts for being inside/outside the sphere
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

ray_color :: proc(r: ray, bounces_remaining: u8, scene: scene) -> float3 {
    if bounces_remaining == 0 do return float3(0.0)

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
        return attenuation * ray_color(new_ray, bounces_remaining - 1, scene)
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

    samples_per_pixel := 50
    max_ray_depth : u8 = 10

    camera := init_camera(
        640,
        360,
        //1280,
        //720,
        float3{13.0, 2.0, 3.0},
        float3{0.0, 0.0, 0.0},
        math.PI / 8.0,
        0.017453292519943295
        //0.0
    )

    //Define scene objects
    scene := scene {
        spheres = [dynamic]sphere{
            {
                origin = float3{0.0, -1000.0, 0.0},
                radius = 1000.0,
                material_idx = 0
            },
            {
                origin = float3{0.0, 1.0, 0.0},
                radius = 1.0,
                material_idx = 1
            },
            {
                origin = float3{-4.0, 1.0, 0.0},
                radius = 1.0,
                material_idx = 2
            },
            {
                origin = float3{4.0, 1.0, 0.0},
                radius = 1.0,
                material_idx = 3
            }
        },
        materials = [dynamic]material {
            {
                albedo = float3(0.5),
                type = .Matte
            },
            {
                refractive_index = 1.5,
                type = .Dielectric
            },
            {
                albedo = float3{0.4, 0.2, 0.1},
                type = .Matte
            },
            {
                albedo = float3{0.7, 0.6, 0.5},
                type = .Metal
            }
        }
    }

    //Randomly add many spheres and materials
    for a := -11; a < 11; a += 1 {
        for b := -11; b < 11; b += 1 {
            rand_mat := rand.float32()
            origin := float3{f32(a) + 0.9*rand.float32(), 0.2, f32(b) + 0.9*rand.float32() }

            if length(origin - float3{4.0, 0.2, 0.0}) > 0.9 {
                sphere := sphere {
                    origin = origin,
                    radius = 0.2,
                    material_idx = u64(len(scene.materials))
                }
                
                mat: material
                if rand_mat < 0.55 {
                    //matte
                    mat.type = .Matte
                    mat.albedo = rand_float3() * rand_float3()
                    append(&scene.spheres, sphere)
                    append(&scene.materials, mat)
                } else if rand_mat < 0.9 {
                    //metal
                    mat.type = .Metal
                    mat.albedo = rand_float3_range(0.5, 1.0)
                    mat.roughness = rand.float32_range(0.0, 0.5)
                    append(&scene.spheres, sphere)
                    append(&scene.materials, mat)
                } else {
                    //glass
                    mat.type = .Dielectric
                    mat.refractive_index = 1.5
                    append(&scene.spheres, sphere)
                    append(&scene.materials, mat)
                }
            }
        }
    }

    //Iterate through the image's pixels
    //from top-left to bottom-right, 
    fmt.println("Beginning main tracing routine...")
    for j := 0; j < camera.image_height; j += 1 {
        fmt.print((camera.image_height - j), "lines remaining...     \r")
        for i := 0; i < camera.image_width; i += 1 {
            //Compute color
            color := float3(0.0)

            //We send multiple, randomly jittered rays through the image plane
            for s := 0; s < samples_per_pixel; s += 1 {
                //Generate the ray to send through the current pixel
                ray := get_pixel_ray(camera, u32(i), u32(j))
    
                //Trace ray against scene
                color += ray_color(ray, max_ray_depth, scene)
            }
            color /= float3(samples_per_pixel)

            //Write color to correct pixel
            pixel_idx := j * camera.image_width + i
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
    fmt.fprintln(image_file, camera.image_width, camera.image_height)
    os.write_string(image_file, "255\n")

    fmt.println("Writing framebuffer to file...")
    for j := 0; j < camera.image_height; j += 1 {
        for i := 0; i < camera.image_width; i += 1 {
            //One can imagine SIMDing the hell out of this
            pixel_idx := j * camera.image_width + i
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