#version 300 es
    precision highp float;

    //time for animation
    uniform float u_time;
    uniform vec2 u_resolution;
    uniform int u_maxBounces;
    uniform float u_ambientStrength;

    // Camera params
    uniform vec3 u_camPos;
    uniform vec3 u_camForward;
    uniform vec3 u_camRight;
    uniform vec3 u_camUp;

    // Element controls
    uniform vec3  u_centerLightColor;
    uniform float u_centerLightBrightness;

    out vec4 fragColor;

    // Material types
    const int MATERIAL_DIFFUSE = 0;
    const int MATERIAL_REFLECTIVE = 1;
    const int MATERIAL_REFRACTIVE = 2;
    const int MATERIAL_LIGHT = 3;

    // Ray structure
    struct Ray {
        vec3 origin;
        vec3 direction;
    };

    // Sphere structure
    struct Sphere {
        vec3 center;
        float radius;
        vec3 color;
        int material;
        float reflectivity;
        float refractiveIndex;
    };


    // Hit record structure
    struct HitData {
        bool hit;
        float t;
        vec3 point;
        vec3 normal;
        vec3 color;
        int material;
        float reflectivity;
        float refractiveIndex;
        bool frontFace;
    };

    // Scene objects
    Sphere spheres[1];          // glass sphere
    Sphere lightSphere;

    // Initialize scene objects
    void initScene(float time) {

        // Center light sphere
        vec3 centerCol = clamp(u_centerLightColor, 0.0, 1.0) * u_centerLightBrightness;
        lightSphere = Sphere(
            vec3(0.0),
            1.5,
            centerCol,
            MATERIAL_LIGHT,
            0.0,
            1.0
        );

        // Glass sphere
        float glassOrbitR = 3.0;
        float glassAngle = time * 0.4;
        vec3 sCenter = vec3(
            glassOrbitR * cos(glassAngle + 2.1),
            0.0,
            glassOrbitR * sin(glassAngle + 2.1)
        );
        spheres[0] = Sphere(
            sCenter,
            0.8,
            vec3(0.9, 0.95, 1.0),
            MATERIAL_REFRACTIVE,
            0.1,
            1.5
        );
    }

    bool intersectSphere(Ray ray, Sphere sphere, out float t) {
        vec3 oc = ray.origin - sphere.center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0 * dot(oc, ray.direction);
        float c = dot(oc, oc) - sphere.radius * sphere.radius;
        float discriminant = b * b - 4.0 * a * c;

        if (discriminant < 0.0) {
            return false;
        }

        float sqrtD = sqrt(discriminant);
        float t0 = (-b - sqrtD) / (2.0 * a);
        float t1 = (-b + sqrtD) / (2.0 * a);

        t = t0;
        if (t < 0.001) {
            t = t1;
            if (t < 0.001) {
                return false;
            }
        }
        return true;
    }

    HitData traceScene(Ray ray, bool includeLightSphere) {
        HitData closest;
        closest.hit = false;
        closest.t = 1e10;

        // Light sphere intersection
        if (includeLightSphere) {
            float t;
            if (intersectSphere(ray, lightSphere, t)) {
                if (t < closest.t) {
                    closest.hit = true;
                    closest.t = t;
                    closest.point = ray.origin + t * ray.direction;
                    vec3 n = normalize(closest.point - lightSphere.center);
                    closest.frontFace = dot(ray.direction, n) < 0.0;
                    closest.normal = closest.frontFace ? n : -n;
                    closest.color = lightSphere.color;
                    closest.material = lightSphere.material;
                    closest.reflectivity = lightSphere.reflectivity;
                    closest.refractiveIndex = lightSphere.refractiveIndex;
                }
            }
        }

        // Glass sphere intersection
        for (int i = 0; i < 1; i++) {
            float t;
            if (intersectSphere(ray, spheres[i], t)) {
                if (t < closest.t) {
                    closest.hit = true;
                    closest.t = t;
                    closest.point = ray.origin + t * ray.direction;
                    vec3 n = normalize(closest.point - spheres[i].center);
                    closest.frontFace = dot(ray.direction, n) < 0.0;
                    closest.normal = closest.frontFace ? n : -n;
                    closest.color = spheres[i].color;
                    closest.material = spheres[i].material;
                    closest.reflectivity = spheres[i].reflectivity;
                    closest.refractiveIndex = spheres[i].refractiveIndex;
                }
            }
        }

        return closest;
    }

    // refraction
    vec3 refraction(vec3 I, vec3 N, float eta) {
        return refract(I, N, eta);
    }

    // Schlick
    float schlick(float cosine, float ref_idx) {
        float r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
    }

    vec3 trace(Ray ray, int maxDepth) {
        vec3 color = vec3(0.0);

        vec3 attenuation = vec3(1.0);

        for (int depth = 0; depth < 8; depth++) {
            if (depth >= maxDepth) break;

            HitData lightHit = traceScene(ray, true);
            // First object hit is a light, accumulate emission then stop
            if (lightHit.hit && lightHit.material == MATERIAL_LIGHT) {
                color += attenuation * lightHit.color;
                break;
            }

            HitData hit = traceScene(ray, false);

            if (!hit.hit) {
                // No hit
                vec3 skyColor;

                // Use sky gradient
                vec3 unitDir = normalize(ray.direction);
                float tEnv = 0.5 * (unitDir.y + 1.0);
                skyColor = mix(vec3(0.05, 0.05, 0.1), vec3(0.4, 0.7, 1.0), tEnv);
               
                color += attenuation * skyColor;
                break;
            }

            // Central light contribution
            vec3 viewDir = normalize(-ray.direction);
            vec3 lightDir1 = normalize(lightSphere.center - hit.point);
            float lightDist1 = length(lightSphere.center - hit.point);

            Ray shadowRay1 = Ray(hit.point + lightDir1 * 0.001, lightDir1);
            HitData shadowHit1 = traceScene(shadowRay1, false);
            bool inShadow1 = shadowHit1.hit && shadowHit1.t < lightDist1;

            float diff1 = max(dot(hit.normal, lightDir1), 0.0);
            vec3 halfDir1 = normalize(lightDir1 + viewDir);
            float spec1 = pow(max(dot(hit.normal, halfDir1), 0.0), 64.0);
            vec3 lightColor1 = lightSphere.color;

            vec3 contrib1 = vec3(0.0);
            if (!inShadow1) {
                contrib1 = hit.color * diff1 * lightColor1 + spec1 * lightColor1;
            }

            vec3 directLight = contrib1;
            vec3 ambientLight = u_ambientStrength * hit.color;

            if (hit.material == MATERIAL_DIFFUSE) {
                // local illumination only
                vec3 lighting = directLight + ambientLight;
                color += attenuation * lighting;
                break;

            } else if (hit.material == MATERIAL_REFRACTIVE) {
                // Refract or reflect
                float etaRatio = hit.frontFace ? (1.0 / hit.refractiveIndex) : hit.refractiveIndex;
                float cosTheta = min(dot(-ray.direction, hit.normal), 1.0);
                float sin2Theta = 1.0 - cosTheta * cosTheta;
                bool total_internel_reflectioin = etaRatio * etaRatio * sin2Theta > 1.0;
                float schlick_value = schlick(cosTheta, hit.refractiveIndex);
                vec3 direction;
                if (total_internel_reflectioin || schlick_value > 0.5) {
                    direction = reflect(ray.direction, hit.normal);
                } else {
                    direction = refraction(ray.direction, hit.normal, etaRatio);
                }

                ray = Ray(hit.point + direction * 0.001, direction);
                attenuation *= hit.color;

            } else {
                // Mirror surface: mostly reflection with a bit of local lighting
                color += attenuation * (0.02 * directLight + 0.02 * ambientLight);

                attenuation *= hit.reflectivity;
                ray = Ray(hit.point + hit.normal * 0.001, reflect(ray.direction, hit.normal));
            }

            if (length(attenuation) < 0.01) break;
        }

        return color;
    }

    void main() {
        initScene(u_time);

        vec2 uv = (gl_FragCoord.xy - 0.5 * u_resolution) / u_resolution.y;

        // Build ray from camera
        float focalDist = 1.5;
        vec3 direction = normalize(
            u_camForward * focalDist +
            u_camRight * uv.x +
            u_camUp * uv.y
        );

        Ray ray = Ray(u_camPos, direction);
        vec3 color = trace(ray, u_maxBounces);

        color = pow(color, vec3(1.0 / 2.2));

        fragColor = vec4(color, 1.0);
    }