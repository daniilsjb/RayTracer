#include <SFML/Graphics.hpp>

#include <math.h>
#include <vector>

struct Vector3
{
    float x, y, z;

    float Magnitude() const
    {
        return sqrtf(x * x + y * y + z * z);
    }

    void Normalize()
    {
        float magnitude = Magnitude();

        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
    }

    Vector3 Normalized() const
    {
        float magnitude = Magnitude();
        return { x / magnitude, y / magnitude, z / magnitude };
    }

    Vector3 operator-() const
    {
        return { -x, -y, -z };
    }

    Vector3 operator+(const Vector3& other) const
    {
        return { x + other.x, y + other.y, z + other.z };
    }

    Vector3 operator-(const Vector3& other) const
    {
        return { x - other.x, y - other.y, z - other.z };
    }

    Vector3 operator*(const Vector3& other) const
    {
        return { x * other.x, y * other.y, z * other.z };
    }

    Vector3 operator/(const Vector3& other) const
    {
        return { x / other.x, y / other.y, z / other.z };
    }

    Vector3 operator*(float scalar) const
    {
        return { x * scalar, y * scalar, z * scalar };
    }

    Vector3 operator/(float scalar) const
    {
        return { x / scalar, y / scalar, z / scalar };
    }

    void operator+=(const Vector3& other)
    {
        *this = *this + other;
    }

    void operator-=(const Vector3& other)
    {
        *this = *this - other;
    }

    void operator*=(float scalar)
    {
        *this = *this * scalar;
    }

    void operator/=(float scalar)
    {
        *this = *this / scalar;
    }
};

float Dot(const Vector3& lhs, const Vector3& rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vector3 Reflect(const Vector3& direction, const Vector3& normal)
{
    return direction - (normal * 2.0f * Dot(direction, normal));
}

//Source: https://en.wikipedia.org/wiki/Snell's_law#Vector_form
Vector3 Refract(const Vector3& direction, const Vector3& normal, const float refractionIndex)
{
    //Find the cosine of the angle of incindence
    //Negate the result since ray and normal point in different directions
    float c = -std::max(-1.0f, std::min(Dot(direction, normal), 1.0f));

    //Define refractive indices (air for the outer medium, surface material for the inner medium)
    float n1 = 1.0f;
    float n2 = refractionIndex;

    Vector3 n = normal;

    //If ray is inside the object, then the outer and inner medium need to be swapped
    if (c < 0.0f)
    {
        c = -c;
        std::swap(n1, n2);
        n = -normal;
    }

    float r = n1 / n2;

    //Find the cosine of the angle of refraction
    float k = 1.0f - r * r * (1.0f - c * c);

    //If the angle was greater than the critical angle, then there is no refraction
    if (k < 0.0f)
        return { 0.0f, 0.0f, 0.0f };

    return direction * r + n * (r * c - sqrtf(k));
}

struct Material
{
    Vector3 ambient;
    Vector3 diffuse;
    Vector3 specular;
    float shininess;
    Vector3 reflection;
    Vector3 refraction;
    float refractiveIndex;
};

struct Ray
{
    Vector3 origin;
    Vector3 direction;
};

struct RayHit
{
    float distance;
    Vector3 point;
    Vector3 normal;
    Material material;
};

struct Object
{
    Object() {}
    virtual ~Object() {}

    virtual bool Intersects(const Ray& ray, RayHit& hit) const = 0;
};

struct Plane : public Object
{
    Vector3 normal;
    float d;
    Material material;

    Plane(const Vector3& normal, const float d, const Material& material) : normal(normal), d(d), material(material) {}
    ~Plane() {}

    //Source: https://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld017.htm
    bool Intersects(const Ray& ray, RayHit& hit) const override
    {
        float denom = Dot(ray.direction, normal);

        //Prevent division-by-zero errors
        if (fabs(denom) <= FLT_EPSILON)
            return false;

        //The numerator is flipped to better align parameter d with coordinate axes
        float t = (d - Dot(ray.origin, normal)) / denom;

        if (t <= FLT_EPSILON)
            return false;

        hit.distance = t;
        hit.normal =  (denom > 0.0f) ? -normal : normal;
        hit.point = ray.origin + ray.direction * hit.distance;
        hit.material = material;

        return true;
    }
};

struct Sphere : public Object
{
    Vector3 origin;
    float radius;
    Material material;

    Sphere(const Vector3& origin, const float radius, const Material& material) : origin(origin), radius(radius), material(material) {}
    ~Sphere() {}

    //Source: http://kylehalladay.com/blog/tutorial/math/2013/12/24/Ray-Sphere-Intersection.html
    bool Intersects(const Ray& ray, RayHit& hit) const override
    {
        Vector3 L = origin - ray.origin;
        float tc = Dot(L, ray.direction);
        float dSqr = Dot(L, L) - tc * tc;

        float rSqr = radius * radius;
        if (dSqr > rSqr)
            return false;

        float th = sqrtf(rSqr - dSqr);

        float t0 = tc - th;
        float t1 = tc + th;

        //We're only interested in the closest intersection point
        float t = t0;

        //However, the point is not a valid intersection if it's negative
        if (t < 0.0f)
            t = t1;

        //If both are negative, no intersection actually occured
        if (t < 0.0f)
            return false;

        hit.distance = t;
        hit.point = ray.origin + (ray.direction * hit.distance);
        hit.normal = (hit.point - origin).Normalized();
        hit.material = material;

        return true;
    }
};

struct LightSource
{
    Vector3 origin;
    Vector3 color;
    float intensity;
};

struct Scene
{
    std::vector<Object*> objects;
    std::vector<LightSource> lightSources;
    Vector3 ambientColor;
    Vector3 backgroundColor;

    Scene() : ambientColor { 1.0f, 1.0f, 1.0f }, backgroundColor { 1.0f, 1.0f, 1.0f } {}
    ~Scene()
    {
        for (const auto& object : objects)
            delete object;
    }

    bool CastRay(const Ray& ray, RayHit& hit)
    {
        //Assume that we've hit an object that is as far away as possible
        hit.distance = FLT_MAX;

        //Find objects that intersect with the ray and choose the one that's closest
        for (const auto& object : objects)
        {
            RayHit currentHit;
            if (object->Intersects(ray, currentHit) && currentHit.distance < hit.distance)
                hit = currentHit;
        }

        //We've hit something only if our initial assumption was incorrect
        return hit.distance < FLT_MAX;
    };

    Vector3 TraceRay(const Ray& ray, int currentDepth, const int& maxDepth)
    {
        //Find the closest intersection
        RayHit hit;
        if (currentDepth >= maxDepth || !CastRay(ray, hit))
            return backgroundColor;

        //Retrieve some data about the ray intersection for convenience
        const Material& material = hit.material;
        const Vector3& normal = hit.normal;
        const Vector3& point = hit.point;

        //Accumulate light intensities for all reflective components for this point
        Vector3 ambientIntensity = Vector3 { 1.0f, 1.0f, 1.0f } * ambientColor;
        Vector3 diffuseIntensity = Vector3 { 0.0f, 0.0f, 0.0f };
        Vector3 specularIntensity = Vector3 { 0.0f, 0.0f, 0.0f };

        //Source: https://en.wikipedia.org/wiki/Phong_reflection_model
        for (const auto& lightSource : lightSources)
        {
            const Vector3 lightDirection = (lightSource.origin - point).Normalized();
            const float lightDistance = (lightSource.origin - point).Magnitude();

            //If another object blocks the path from this light source, then it shouldn't be applied
            RayHit obstacleHit;
            if (CastRay({ point + lightDirection * 0.01f, lightDirection }, obstacleHit) && obstacleHit.distance < lightDistance)
                continue;

            const Vector3 reflection = Reflect(lightDirection, normal).Normalized();

            diffuseIntensity += lightSource.color * lightSource.intensity * std::max(Dot(lightDirection, normal), 0.0f);
            specularIntensity += lightSource.color * lightSource.intensity * std::powf(std::max(Dot(reflection, ray.direction), 0.0f), material.shininess);
        }

        //Combine material information with accumulated intensities
        const Vector3 ambient = material.ambient * ambientIntensity;
        const Vector3 diffuse = material.diffuse * diffuseIntensity;
        const Vector3 specular = material.specular * specularIntensity;

        //Combine resulting lighting colors
        Vector3 pointColor = ambient + diffuse + specular;

        //Perform recursive reflection
        const Vector3 reflectedDirection = Reflect(ray.direction, normal).Normalized();
        const Ray reflectedRay = { point + reflectedDirection * 0.01f, reflectedDirection };
        Vector3 reflectColor = material.reflection * TraceRay(reflectedRay, currentDepth + 1, maxDepth);

        //Perform recursive refraction
        const Vector3 refractedDirection = Refract(ray.direction, normal, material.refractiveIndex);
        const Ray refractedRay = { point + refractedDirection * 0.01f, refractedDirection };
        Vector3 refractColor = material.refraction * TraceRay(refractedRay, currentDepth + 1, maxDepth);

        //Combine colors and clamp them to valid range
        Vector3 combinedColor = (pointColor + reflectColor + refractColor);
        combinedColor.x = std::min(combinedColor.x, 1.0f);
        combinedColor.y = std::min(combinedColor.y, 1.0f);
        combinedColor.z = std::min(combinedColor.z, 1.0f);

        return combinedColor;
    }
};

void Render(const int width, const int height, const Vector3* colors);

int main()
{
    //Define screen
    const int screenWidth = 1200;
    const int screenHeight = 800;
    const float aspectRatio = (float)screenWidth / (float)screenHeight;
    Vector3* colors = new Vector3[screenWidth * screenHeight];

    //Defube camera
    const Vector3 cameraOrigin = { 0.0f, 0.0f, -1.0f };
    const float FoV = 3.142592f / 2.0f;
    const int maxDepth = 4;

    //Taken from here: http://devernay.free.fr/cours/opengl/materials.html
    Material emerald;             emerald.ambient = { 0.0215f, 0.1745f, 0.0215f };          emerald.diffuse = { 0.07568f, 0.61424f, 0.07568f };         emerald.specular = { 0.633f, 0.727811f, 0.633f };                   emerald.shininess = 76.8f;         emerald.reflection = { 0.2f, 0.2f, 0.2f };          emerald.refraction = { 0.0f, 0.0f, 0.0f };       emerald.refractiveIndex = 1.0f;
    Material redRubber;         redRubber.ambient = { 0.05f, 0.0f, 0.0f };                redRubber.diffuse = { 0.5f, 0.4f, 0.4f };                   redRubber.specular = { 0.7f, 0.04f, 0.04f };                        redRubber.shininess = 10.0f;       redRubber.reflection = { 0.0f, 0.0f, 0.0f };        redRubber.refraction = { 0.0f, 0.0f, 0.0f };     redRubber.refractiveIndex = 1.0f;
    Material bronze;               bronze.ambient = { 0.2125f, 0.1275f, 0.054f };            bronze.diffuse = { 0.714f, 0.4284f, 0.18144f };             bronze.specular = { 0.393548f, 0.271906f, 0.166721f };              bronze.shininess = 25.6f;          bronze.reflection = { 0.05f, 0.05f, 0.05f };        bronze.refraction = { 0.0f, 0.0f, 0.0f };        bronze.refractiveIndex = 1.0f;
    Material pearl;                 pearl.ambient = { 0.25f, 0.20725f, 0.20725f };            pearl.diffuse = { 1.0f, 0.829f, 0.829f };                   pearl.specular = { 0.296648f, 0.296648f, 0.296648f };               pearl.shininess = 225.28f;         pearl.reflection = { 0.1f, 0.1f, 0.1f };            pearl.refraction = { 0.0f, 0.0f, 0.0f };         pearl.refractiveIndex = 1.0f;
    Material ruby;                   ruby.ambient = { 0.1745f, 0.01175f, 0.01175f };           ruby.diffuse = { 0.61424f, 0.04136f, 0.04136f };            ruby.specular = { 0.727811f, 0.626959f, 0.626959f };                ruby.shininess = 76.8f;            ruby.reflection = { 0.2f, 0.2f, 0.2f };             ruby.refraction = { 0.0f, 0.0f, 0.0f };          ruby.refractiveIndex = 1.0f;
    Material cyanPlastic;     cyanPlastic.ambient = { 0.0f, 0.1f, 0.06f };              cyanPlastic.diffuse = { 0.0f, 0.50980392f, 0.50980392f };   cyanPlastic.specular = { 0.50196078f, 0.50196078f, 0.50196078f };   cyanPlastic.shininess = 32.0f;     cyanPlastic.reflection = { 0.03f, 0.03f, 0.03f };   cyanPlastic.refraction = { 0.0f, 0.0f, 0.0f };   cyanPlastic.refractiveIndex = 1.0f;
    Material yellowPlastic; yellowPlastic.ambient = { 0.0f, 0.0f, 0.0f };             yellowPlastic.diffuse = { 0.6f, 0.6f, 0.5f };               yellowPlastic.specular = { 0.6f, 0.6f, 0.5f };                      yellowPlastic.shininess = 32.0f;   yellowPlastic.reflection = { 0.03f, 0.03f, 0.03f }; yellowPlastic.refraction = { 0.0f, 0.0f, 0.0f }; yellowPlastic.refractiveIndex = 1.0f;
    Material redPlastic;       redPlastic.ambient = { 0.0f, 0.0f, 0.0f };                redPlastic.diffuse = { 0.5f, 0.0f, 0.0f };                  redPlastic.specular = { 0.7f, 0.6f, 0.6f };                         redPlastic.shininess = 32.0f;      redPlastic.reflection = { 0.03f, 0.03f, 0.03f };    redPlastic.refraction = { 0.0f, 0.0f, 0.0f };    redPlastic.refractiveIndex = 1.0f;
    Material silver;               silver.ambient = { 0.19225f, 0.19225f, 0.19225f };        silver.diffuse = { 0.50754f, 0.50754f, 0.50754f };          silver.specular = { 0.508273f, 0.508273f, 0.508273f };              silver.shininess = 51.2f;          silver.reflection = { 0.2f, 0.2f, 0.2f };           silver.refraction = { 0.0f, 0.0f, 0.0f };        silver.refractiveIndex = 1.0f;
    Material mirror;               mirror.ambient = { 0.0f, 0.0f, 0.0f };                    mirror.diffuse = { 0.0f, 0.0f, 0.0f };                      mirror.specular = { 1.0f, 1.0f, 1.0f };                             mirror.shininess = 1500.0f;        mirror.reflection = { 0.8f, 0.8f, 0.8f };           mirror.refraction = { 0.0f, 0.0f, 0.0f };        mirror.refractiveIndex = 1.0f;
    Material glass;                 glass.ambient = { 0.0f, 0.0f, 0.0f };                     glass.diffuse = { 0.0f, 0.0f, 0.0f };                       glass.specular = { 0.3f, 0.3f, 0.3f };                              glass.shininess = 75.0f;           glass.reflection = { 0.1f, 0.1f, 0.1f };            glass.refraction = { 0.8f, 0.8f, 0.8f };         glass.refractiveIndex = 1.5f;

    //Describe the scene that will be rendered
    Scene scene;
    scene.lightSources = {
        { { 15.0f, 30.0f, 13.0f }, { 1.0f, 1.0f, 1.0f }, 1.0f },
        { { 8.0f, 20.0f, -5.0f }, { 0.8f, 0.7f, 0.6f }, 0.7f },
        { { -14.0f, 9.0f, 12.0f }, { 1.0f, 1.0f, 1.0f}, 0.3f },
    };
    scene.objects = {
        new Sphere({ 10.0f, 0.0f, 30.0f }, 6.0f, emerald),
        new Sphere({ -12.0f, -2.0f, 15.0f }, 4.0f, redRubber),
        new Sphere({ 20.0f, 20.0f, 35.0f }, 10.0f, bronze),
        new Sphere({ -3.0f, 3.0f, 18.0f }, 5.0f, pearl),
        new Sphere({ 0.0f, -2.0f, 20.0f }, 4.0f, mirror),
        new Sphere({ -15.0f, 8.0f, 16.0f }, 2.0f, mirror),
        new Sphere({ 5.0f, 10.0f, 25.0f }, 3.0f, glass),
        new Sphere({ -3.0f, 0.0f, 4.0f }, 1.0f, glass),
        new Plane({ 0.0f, 1.0f, 0.0f }, -15.0f, silver)
    };
    scene.backgroundColor = { 0.2f, 0.8f, 0.9f };

    //Perform ray tracing for each pixel
    for (int y = 0; y < screenHeight; y++)
    {
        for (int x = 0; x < screenWidth; x++)
        {
            //Normalize the center of each pixel to range [-1; 1]
            const float normalizedX = 2.0f * (x + 0.5f) / ((float)screenWidth  - 1.0f) - 1.0f;
            const float normalizedY = 2.0f * (y + 0.5f) / ((float)screenHeight - 1.0f) - 1.0f;

            //Map the normalized coordinates to a direction within Field of View
            Vector3 direction = {
                normalizedX * std::tan(FoV * 0.5f) * aspectRatio,
               -normalizedY * std::tan(FoV * 0.5f),
                1.0f
            };

            direction.Normalize();

            //Trace a ray from the camera in computed direction
            const Ray ray = { cameraOrigin, direction };
            colors[y * screenWidth + x] = scene.TraceRay(ray, 0, maxDepth);
        }
    }

    Render(screenWidth, screenHeight, colors);

    delete[] colors;
}

void Render(const int width, const int height, const Vector3* colors)
{
    sf::RenderWindow window(sf::VideoMode(width, height), "Ray Tracing");

    sf::Image targetImage;
    targetImage.create(width, height);
    for (unsigned int y = 0; y < (unsigned int)height; y++)
    {
        for (unsigned int x = 0; x < (unsigned int)width; x++)
        {
            const Vector3& color = colors[y * width + x];
            targetImage.setPixel(x, y, sf::Color((sf::Uint8)(color.x * 255.0f), (sf::Uint8)(color.y * 255.0f), (sf::Uint8)(color.z * 255.0f)));
        }
    }

    targetImage.saveToFile("output.png");

    sf::Texture targetTexture;
    targetTexture.loadFromImage(targetImage);

    sf::Sprite targetSprite;
    targetSprite.setTexture(targetTexture);

    while (window.isOpen())
    {
        sf::Event ev;
        while (window.pollEvent(ev))
        {
            if (ev.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
        window.draw(targetSprite);
        window.display();
    }
}
