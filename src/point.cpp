/*
Copyright (c) 2022 King's College London, MeTrICS Lab, Renato Besenczi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include "point.h"

namespace newresampler {

void Point::normalize() {
    double n = norm();
    if (n != 0)
    {
        X /= n;
        Y /= n;
        Z /= n;
    }
}

bool SameSide(const Point &p1, const Point &p2, const Point &a, const Point &b) {

    Point cp1 = (b - a) * (p1 - a);
    Point cp2 = (b - a) * (p2 - a);

    if ((cp1 | cp2) >= -1E-6) return true; else return false;
}

bool PointInTriangle(const Point &p, const Point &a, const Point &b, const Point &c) {

    if (SameSide(p, a, b, c) && SameSide(p, b, c, a) && SameSide(p, c, a, b))
        return true;
    else return false;
}

Point projectPoint(const Point &vb, const Point &v1, const Point &v2,
                   const Point &v3, Point &PP) {

    Point s1 = v3 - v1;
    s1.normalize();
    Point s2 = v2 - v1;
    s2.normalize();
    Point s3 = s1 * s2; // s3= normal
    s3.normalize();

    double si = (s3 | v1) / (s3 | vb);
    // formula for line plane intersection s3.(v1-sivb)=0
    // i.e. (v1-sivb) should be perpendicular to plane normal s3

    PP = vb * si;  // PP therefore where vb intersects plane
    return s3;
}

double computeArea(const Point &v0, const Point &v1, const Point &v2) {

    Point res, v0v1, v0v2;

    v0v1 = v1 - v0;
    v0v2 = v2 - v0;
    res = v0v1 * v0v2;

    return 0.5 * res.norm();
}

double operator|(const Point &v1, const Point &v2) {
    return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
}

Point operator*(const Point &v1, const Point &v2) {

    return {v1.Y * v2.Z - v1.Z * v2.Y,
            v2.X * v1.Z - v2.Z * v1.X,
            v1.X * v2.Y - v2.X * v1.Y};
}

Point operator+(const Point &v1, const Point &v2) {
    return {v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z};
}

Point operator-(const Point &v1, const Point &v2) {
    return {v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z};
}

Point operator/(const Point &v, const double &d) {
    if (d != 0)
        return {v.X / d, v.Y / d, v.Z / d};
    else
        throw MeshException{"newresampler::Point operator/ division by zero"};
}

Point operator*(const Point &v, const double &d) {
    return {v.X * d, v.Y * d, v.Z * d};
}

Point operator*(const NEWMAT::Matrix &M, const Point &v) {
    if (M.Ncols() != 3 || M.Nrows() != 3)
        throw MeshException{"newresampler::Point matrix multiply error: matrix should be 3x3"};
    return {M(1, 1) * v.X + M(1, 2) * v.Y + M(1, 3) * v.Z,
            M(2, 1) * v.X + M(2, 2) * v.Y + M(2, 3) * v.Z,
            M(3, 1) * v.X + M(3, 2) * v.Y + M(3, 3) * v.Z};
}

Point operator*(const Point &v, const NEWMAT::Matrix &M) {
    if (M.Ncols() != 3 || M.Nrows() != 3)
        throw MeshException{"newresampler::Point matrix multiply error: matrix should be 3x3"};
    return {v.X * M(1, 1) + v.Y * M(2, 1) + v.Z * M(3, 1),
            v.X * M(1, 2) + v.Y * M(2, 2) + v.Z * M(3, 2),
            v.X * M(1, 3) + v.Y * M(2, 3) + v.Z * M(3, 3)};
}

Point operator+=(Point &p1, const Point &p2) {
    p1.X += p2.X;
    p1.Y += p2.Y;
    p1.Z += p2.Z;
    return p1;
}

Point operator*=(Point &p, const double &d) {
    p.X *= d;
    p.Y *= d;
    p.Z *= d;
    return p;
}

Point operator/=(Point &p, const double &d) {
    if (d != 0) {
        p.X /= d;
        p.Y /= d;
        p.Z /= d;
    } else throw MeshException{"newresampler::Point operator/= division by zero"};
    return p;
}

bool operator==(const Point &p1, const Point &p2) {
    return ((std::fabs(p1.X - p2.X) < 1e-4)
        && (std::fabs(p1.Y - p2.Y) < 1e-4)
        && (std::fabs(p1.Z - p2.Z) < 1e-4));
}

bool operator!=(const Point& p1, const Point& p2) {
    return !(p1==p2);
}

std::ostream& operator<<(std::ostream& os, const Point& pt) {
    return os << pt.X << ' ' << pt.Y << ' ' << pt.Z;
}

} //namespace newresampler