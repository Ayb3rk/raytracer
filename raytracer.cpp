#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>

typedef unsigned char RGB[3];

struct Ray {
    parser::Vec3f origin;
    parser::Vec3f direction;
};

parser::Vec3f multS(parser::Vec3f a,float s)
{
    parser::Vec3f result{};
    result.x = a.x*s;
    result.y = a.y*s;
    result.z = a.z*s;
    return result;
}

float dot(parser::Vec3f a,parser::Vec3f b)
{
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

parser::Vec3f normalize(parser::Vec3f a)
{
    return multS(a,1.0/sqrt(dot(a,a)));
}

parser::Vec3f add(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result;
    result.x = a.x+b.x;
    result.y = a.y+b.y;
    result.z = a.z+b.z;
    return result;
}
parser::Vec3f CrossProduct(parser::Vec3f a, parser::Vec3f b) {
    parser::Vec3f result{};
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

Ray GenerateRay(int i, int j, parser::Camera camera, int width, int height)
{
    Ray result{};
    float su,sv;
    parser::Vec3f m{},q{},s{};

    auto left = camera.near_plane.x;
    auto right = camera.near_plane.y;
    auto bottom = camera.near_plane.z;
    auto top = camera.near_plane.w;

    su = (i+0.5)*(right - left)/width;
    sv = (j+0.5)*(top - bottom)/height;
    parser::Vec3f v = camera.up;
    parser::Vec3f w = multS(camera.gaze,-1);
    parser::Vec3f u = CrossProduct(v,w);

    m = add(camera.position,multS(camera.gaze,camera.near_distance));
    q = add(m, add(multS(u,left),multS(v,top)));
    s = add(q,add(multS(u,su),multS(v,-sv)));
    result.origin = camera.position;
    result.direction = add(s,multS(camera.position,-1));
    return result;
}



double intersectSphere (parser::Vec3f sphereCenter, double sphereRadius, Ray ray) {
    float A,B,C; //constants for the quadratic equation

    float delta;

    parser::Vec3f c;

    c = sphereCenter;

    float t,t1,t2;

    C = (ray.origin.x-c.x)*(ray.origin.x-c.x)+(ray.origin.y-c.y)*(ray.origin.y-c.y)+(ray.origin.z-c.z)*(ray.origin.z-c.z)-sphereRadius*sphereRadius;

    B = 2*ray.direction.x*(ray.origin.x-c.x)+2*ray.direction.y*(ray.origin.y-c.y)+2*ray.direction.z*(ray.origin.z-c.z);

    A = ray.direction.x*ray.direction.x+ray.direction.y*ray.direction.y+ray.direction.z*ray.direction.z;

    delta = B*B-4*A*C;

    if (delta<0) return -1;
    else if (delta==0)
    {
        t = -B / (2*A);
    }
    else
    {
        delta = sqrt(delta);
        A = 2*A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;

        if (t1<t2) t=t1; else t=t2;
    }

    return t;
}

parser::Vec3f ComputeColor(Ray ray, parser::Scene scene)
{
    double t = 1000000;
    parser::Vec3f color{};
    parser::Sphere intersectedSphere{};
    for (auto sphere : scene.spheres) {
        auto sphereCenter = scene.vertex_data[sphere.center_vertex_id-1];
        double temp = intersectSphere(sphereCenter, sphere.radius, ray);
        if (temp > 0 && temp < t) {
            t = temp;
            intersectedSphere = sphere;
        }
    }
    if (t == 1000000) {
        return color;
    }
    auto P = add(ray.origin, multS(ray.direction, t));
    auto L = normalize(add(scene.point_lights[0].position, multS(P, -1)));
    auto N = normalize(add(P, multS(scene.vertex_data[intersectedSphere.center_vertex_id-1], -1)));
    color = multS(scene.materials[intersectedSphere.material_id-1].diffuse, dot(L, N));
    return color;

}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    auto spheres = scene.spheres;

    int width = scene.cameras[0].image_width;
    int height = scene.cameras[0].image_height;

    auto backroundColor = scene.background_color;

    auto* image = new unsigned char [width * height * 3]; // 3 channels

    for (int i = 0; i < width * height; ++i) //initialize image to backround color
    {
        image[i*3] = backroundColor.x;
        image[i*3+1] = backroundColor.y;
        image[i*3+2] = backroundColor.z;
    }

    for (int j = 0; j < height; j++) { //main loop
        for (int i = 0; i < width; i++) {
            Ray ray{};
            ray = GenerateRay(i, j, scene.cameras[0], width, height);
            auto color = ComputeColor(ray, scene);
            image[(j * width + i) * 3] = color.x * 255 + scene.ambient_light.x;
            image[(j * width + i) * 3 + 1] = color.y * 255 + scene.ambient_light.y;
            image[(j * width + i) * 3 + 2] = color.z * 255 + scene.ambient_light.z;
        }
    }

    write_ppm("test.ppm", image, width, height);

    delete[] image;

    return 0;
}
