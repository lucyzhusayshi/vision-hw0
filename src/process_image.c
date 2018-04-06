#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

// CHW format

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in

    // use clamp padding strategy if coords out of bounds
    if (x < 0) x = 0;
    if (x >= im.w) x = im.w - 1;
        
    if (y < 0) y = 0;
    if (y >= im.h) y = im.h - 1;

    if (c < 0) c = 0;
    if (c >= im.c) c = im.c - 1;

    int index = x + y*im.w + c*im.w*im.h;
    return im.data[index];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    if ((x >= 0 && x < im.w) && (y >= 0 && y < im.h) && (c >= 0 && c < im.c)) {
        int index = x + y*im.w + c*im.w*im.h;
        im.data[index] = v;
    }
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    
    memcpy(copy.data, im.data, im.h*im.w*im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    
    for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
            int index = x + y*im.w;
            gray.data[index] = 0.299*im.data[index+ 0*im.w*im.h] + 0.587*im.data[index + 1*im.w*im.h] + 0.114*im.data[index+ 2*im.w*im.h];
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    if (c >= 0 && c < im.c) {
        for (int y = 0; y < im.h; y++) {
            for (int x = 0; x < im.w; x++) {
                im.data[x + y*im.w + c*im.w*im.h] += v;
            }
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.w*im.h*im.c; i++) {
        if (im.data[i] < 0) im.data[i] = 0;
        if (im.data[i] > 1) im.data[i] = 1;
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
            int index = x + y*im.w;
            float r = im.data[index + 0*im.w*im.h];
            float g = im.data[index + 1*im.w*im.h];
            float b = im.data[index + 2*im.w*im.h];

            float min = three_way_min(r, g, b);
            float max = three_way_max(r, g, b);
            float c = max - min;

            float v = max;
            float s = (v != 0) ? c / v : 0;
            float h;
            if (c == 0) { h = 0; } 
            else {
                if (v == r) { h = (g-b) / c; }
                else if (v == g) { h = ((b-r) / c) + 2; }
                else /* v == b */ { h = ((r-g) / c) + 4; }
            }
            h /= 6;
            if (h < 0) h++;

            im.data[index + 0*im.w*im.h] = h;
            im.data[index + 1*im.w*im.h] = s;
            im.data[index + 2*im.w*im.h] = v;
        }
    }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
            int index = x + y*im.w;
            float h = im.data[index + 0*im.w*im.h];
            float s = im.data[index + 1*im.w*im.h];
            float v = im.data[index + 2*im.w*im.h];

            float c = s * v;
            float max = v;
            float min = (v == c) ? 0.0 : v - c;

            float r, g, b;
            
            float h_ = h * 6;
            if (c == 0) // all values are the same
                { r = v; g = v; b = v; } 
            else if (h_ > 5 && h_ < 6) 
                { r = max; g = min; b = ((h_ / 6 - 1) * 6 * c - g) * -1; }
            else if (h_ == 5) /* max = b = r and min = g */ 
                { r = max; b = max; g = min; }
            else if (h_ < 5 && h_ > 4) /**/ 
                { b = max; g = min; r = (h_ - 4) * c + g; }
            else if (h_ == 4) /* max = b and min = r = g */ 
                { b = max; r = min; g = min; }
            else if (h_ < 4 && h_ > 3) /**/ 
                {b = max; r = min; g = ((h_ - 4) * c - r) * -1; }
            else if (h_ == 3) /* max = b = g and min = r */ 
                { b = max; g = max; r = min;}
            else if (h_ < 3 && h_ > 2) /**/ 
                { g = max; r = min; b = (h_ - 2) * c + r; }
            else if (h_ == 2) /**/ 
                { g = max; r = min; b = min; }
            else if (h_ < 2 && h_ > 1) /**/ 
                { g = max; b = min; r = ((h_ - 2) * c - b) * -1; }
            else if (h_ == 1) /**/ 
                { r = max; g = max; b = min; }
            else if (h_ < 1 && h_ > 0) /**/ 
                { r = max; b = min; g = h_ * c + b; }
            else /* if (h_ == 0) */ 
                { r = max; g = min; b = min; }

            im.data[index + 0*im.w*im.h] = r;
            im.data[index + 1*im.w*im.h] = g;
            im.data[index + 2*im.w*im.h] = b;

            // if (r < 0 || r > 1 || b < 0 || b > 1 || g < 0 || g > 1) {
            //     printf("(%d, %d): r %f, b %f, g %f", x, y, r, b, g);
            //     printf(" h %f, s %f, v %f, h_ %f, c %f, min %f, max %f \n", h, s, v, h_, c, min, max);
            // }
        }
    }
}

void scale_image(image im, int c, float v)
{
    if (c >= 0 && c < im.c) {
        for (int y = 0; y < im.h; y++) {
            for (int x = 0; x < im.w; x++) {
                im.data[x + y*im.w + c*im.w*im.h] *= v;
            }
        }
    }
}

void rgb_to_ciexyz(image im) 
{
    for (int pixel_y = 0; pixel_y < im.h; pixel_y++) {
        for (int pixel_x = 0; pixel_x < im.w; pixel_x++) {
            int c1 = pixel_x + pixel_y*im.w;
            int c2 = c1 + im.w*im.h;
            int c3 = c1 + 2*im.w*im.h;

            // gamma decompression
            float r = (im.data[c1] > 0.04045) ? powf(((im.data[c1] + 0.055) / 1.055), 2.4) : im.data[c1] / 12.92;
            float g = (im.data[c2] > 0.04045) ? powf(((im.data[c2] + 0.055) / 1.055), 2.4) : im.data[c2] / 12.92;
            float b = (im.data[c3] > 0.04045) ? powf(((im.data[c3] + 0.055) / 1.055), 2.4) : im.data[c3] / 12.92;

            

            float x = 0.4124*r + 0.3576*g + 0.1805*b;
            float y = 0.2126*r + 0.7152*g + 0.0722*b;
            float z = 0.0193*r + 0.1192*g + 0.9505*b;

            im.data[c1] = x;
            im.data[c2] = y;
            im.data[c3] = z;
        }
    }
}

void ciexyz_to_rgb(image im) 
{
    for (int pixel_y = 0; pixel_y < im.h; pixel_y++) {
        for (int pixel_x = 0; pixel_x < im.w; pixel_x++) {
            int c1 = pixel_x + pixel_y*im.w;
            int c2 = c1 + im.w*im.h;
            int c3 = c1 + 2*im.w*im.h;

            float x = im.data[c1];
            float y = im.data[c2];
            float z = im.data[c3];

            float r = 3.2406*x + (-1.5372)*y + (-0.4986)*z;
            float g = (-0.9689)*x + 1.8758*y + 0.0415*z;
            float b = 0.0557*x + (-0.2040)*y + 1.0570*z;

            // gamma compression
            r = (r > 0.0031308) ? 1.055 * powf(r, (float)1.0/2.4) - 0.055 : 12.92 * r;
            g = (g > 0.0031308) ? 1.055 * powf(g, (float)1.0/2.4) - 0.055 : 12.92 * g;
            b = (b > 0.0031308) ? 1.055 * powf(b, (float)1.0/2.4) - 0.055 : 12.92 * b;

            im.data[c1] = r;
            im.data[c2] = g;
            im.data[c3] = b;
        }
    }
}

void ciexyz_to_cieluv(image im) {
    // using white reference from https://en.wikipedia.org/wiki/Illuminant_D65
    float x_n = 95.047;
    float y_n = 100.00;
    float z_n = 108.883;

    float u__n = 4*x_n / (x_n + 15*y_n + 3*z_n);
    float v__n = 9*y_n / (x_n + 15*y_n + 3*z_n);

    for (int pixel_y = 0; pixel_y < im.h; pixel_y++) {
        for (int pixel_x = 0; pixel_x < im.w; pixel_x++) {
            int c1 = pixel_x + pixel_y*im.w;
            int c2 = c1 + im.w*im.h;
            int c3 = c1 + 2*im.w*im.h;

            float x = im.data[c1];
            float y = im.data[c2];
            float z = im.data[c3];

            float denom = x + 15*y + 3*z;
            float u_ = (denom == 0) ? 0 : 4*x / (x + 15*y + 3*z);
            float v_ = (denom == 0) ? 0 : 9*y / (x + 15*y + 3*z);

            float l = (y/y_n > powf((float)6.0/29, 3)) ? 116 * powf(y/y_n, (float)1/3) - 16: powf((float)29.0/3, 3) * y/y_n;
            float u = 13*l*(u_ - u__n);
            float v = 13*l*(v_ - v__n);

            im.data[c1] = l;
            im.data[c2] = u;
            im.data[c3] = v;
        }
    }
}

void cieluv_to_ciexyz(image im) {
    // using white reference from https://en.wikipedia.org/wiki/Illuminant_D65
    float x_n = 95.047;
    float y_n = 100.00;
    float z_n = 108.883;

    float u__n = 4*x_n / (x_n + 15*y_n + 3*z_n);
    float v__n = 9*y_n / (x_n + 15*y_n + 3*z_n);

    for (int pixel_y = 0; pixel_y < im.h; pixel_y++) {
        for (int pixel_x = 0; pixel_x < im.w; pixel_x++) {
            int c1 = pixel_x + pixel_y*im.w;
            int c2 = c1 + im.w*im.h;
            int c3 = c1 + 2*im.w*im.h;

            float l = im.data[c1];
            float u = im.data[c2];
            float v = im.data[c3];

            float u_ = u / (13*l) + u__n;
            float v_ = v / (13*l) + v__n;

            float y = (l > 8) ? y_n * powf((l + 16) / 116, 3) : y_n * l * powf((float)3.0/29, 3);
            float x = y * 9*u_ / (4*v_);
            float z = y * (12 - 3*u_ - 20*v_) / (4*v_);

            im.data[c1] = x;
            im.data[c2] = y;
            im.data[c3] = z;
        }
    }
}

void cieluv_to_hcl(image im) {
    for (int pixel_y = 0; pixel_y < im.h; pixel_y++) {
        for (int pixel_x = 0; pixel_x < im.w; pixel_x++) {
            int c1 = pixel_x + pixel_y*im.w;
            int c2 = c1 + im.w*im.h;
            int c3 = c1 + 2*im.w*im.h;

            float l = im.data[c1];
            float u = im.data[c2];
            float v = im.data[c3];

            float c = sqrtf(powf(u, 2) + powf(v, 2));
            float h = (c != 0) ? atan2f(v, u) : 0;
            
            im.data[c1] = h;
            im.data[c2] = c;
            im.data[c3] = l;
        }
    }
}

void hcl_to_cieluv(image im) {
    for (int pixel_y = 0; pixel_y < im.h; pixel_y++) {
        for (int pixel_x = 0; pixel_x < im.w; pixel_x++) {
            int c1 = pixel_x + pixel_y*im.w;
            int c2 = c1 + im.w*im.h;
            int c3 = c1 + 2*im.w*im.h;

            float h = im.data[c1];
            float c = im.data[c2];
            float l = im.data[c3];

            float u, v;
            if (c == 0) {
                u = 0;
                v = 0;
            } else {
                float theta = tanf(h);
                u = c * sinf(theta);
                v = c * cosf(theta);
            }

            im.data[c1] = l;
            im.data[c2] = u;
            im.data[c3] = v;
        }
    }
}

void rgb_to_hcl(image im) {
    rgb_to_ciexyz(im);
    ciexyz_to_cieluv(im);
    cieluv_to_hcl(im);
}

void hcl_to_rgb(image im) {
    hcl_to_cieluv(im);
    cieluv_to_ciexyz(im);
    ciexyz_to_rgb(im);
    // clamp_image(im);
}