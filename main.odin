package main

import "core:fmt"
import "core:os"
import "vendor:stb/image"

float3 :: [3]f32

ray :: struct {
    origin : float3,
    direction : float3
}

main :: proc() {
    fmt.println("Initiating swag mode...");
    ctxt := context

    image_width := 640
    image_height := 480
    aspect_ratio := f32(image_width) / f32(image_height)

    viewport_height : f32 = 1.0
    viewport_width := aspect_ratio * viewport_height

    //Internal framebuffer we hold in RAM during rendering
    framebuffer : [dynamic]float3
    resize(&framebuffer, image_height * image_width)

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

    for j := 0; j < image_height; j += 1 {
        fmt.print((image_height - j),"lines remaining...     \r")
        for i := 0; i < image_width; i += 1 {
            red := f32(i) / (f32(image_width) - 1.0)
            green := f32(j) / (f32(image_height) - 1.0)
            red_int := i32(red * 255.0)
            green_int := i32(green * 255.0)
            fmt.fprintln(image_file, red_int, green_int, "0", flush = false)


        }
    }
    fmt.println()

    fmt.println("Done!")
}