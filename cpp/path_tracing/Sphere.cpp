#include "Sphere.h"

double Sphere::Intersect(const Ray &ray) const {
    /* Check for ray-sphere intersection by solving for t:
        t^2*d.d + 2*t*(o-p).d + (o-p).(o-p) - R^2 = 0 */
    Vector op = position - ray.org;
    double eps = 1e-4;
    double b = op.dotProduct(ray.dir);
    double radicant = b*b - op.dotProduct(op) + radius*radius;
    if (radicant < 0.0)
        return 0.0;      /* No intersection */
    else
        radicant = sqrt(radicant);

    double t;
    t = b - radicant;    /* Check smaller root first */
    if(t > eps)
        return t;

    t = b + radicant;
    if(t > eps)          /* Check second root */
        return t;

    return 0.0;          /* No intersection in ray direction */
}
