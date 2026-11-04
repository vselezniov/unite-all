#ifndef _COMPARE_FLOAT_H
#define	_COMPARE_FLOAT_H

// Максимальная погрешность при сравнении |double|
extern const double EPSILON;
extern const double SQR_EPSILON;    // квадрат |EPSILON| для сравнения площадей [?]

// Функции сравнения, учитывающие погрешность
inline int compare_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    if (x + eps < y)
        return -1;
    else if (x > y + eps)
        return 1;
    else
        return 0;
}

// Равно
inline bool equal_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return ( (x - eps <= y) && (y <= x + eps) );
}

// Больше
inline bool greater_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x > y + eps;
}

// Больше или равно
inline bool greater_eq_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x >= y - eps;
}

// Меньше
inline bool less_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x < y - eps;
}

// Меньше или равно
inline bool less_eq_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x <= y + eps;
}

// Сравнение точек на плоскости
inline int compare_eps(const leda_point& a, const leda_point& b, const double &eps)
{
    int r = compare_eps(a.xcoord(), b.xcoord(), eps);
    return ( r != 0 ? r : compare_eps(a.ycoord(), b.ycoord(), eps) );
}

inline int compare_eps(const leda_point& a, const leda_point& b)
{
    return compare_eps(a, b, EPSILON);
}

inline int equal_eps(const leda_point& a, const leda_point& b, const double &eps = EPSILON)
{
    return equal_eps(a.xcoord(), b.xcoord(), eps) && equal_eps(a.ycoord(), b.ycoord(), eps);
}

// Вычисляет ориентацию трех точек, учитывая погрешность
inline int orientation_eps(const leda_point& a, const leda_point& b, const leda_point& c,
        const double &eps = SQR_EPSILON)
{
    double d1 = (a.xcoord() - b.xcoord()) * (a.ycoord() - c.ycoord());
    double d2 = (a.ycoord() - b.ycoord()) * (a.xcoord() - c.xcoord());
    return compare_eps(d1, d2, eps);
}

// Вычисляет ориентацию точки относительно направленного отрезка, учитывая погрешность
inline int orientation_eps(const leda_segment& s, const leda_point& p, const double &eps = SQR_EPSILON)
{
    return compare_eps( s.dy()*(s.xcoord1()-p.xcoord()), s.dx()*(s.ycoord1()-p.ycoord()), eps );
}

// Лежат ли три точки на одной прямой
inline bool collinear_eps(const leda_point& a, const leda_point& b, const leda_point& c,
        const double &eps = SQR_EPSILON)
{
    return equal_eps( (a.ycoord()-b.ycoord()) * (a.xcoord()-c.xcoord()),
            (a.xcoord()-b.xcoord()) * (a.ycoord()-c.ycoord()), eps );
}

// Содержит ли сегмент точку
inline bool contains_eps(const leda_segment &s, const leda_point& p)
{
    leda_point a = s.source();
    leda_point b = s.target();
    if(equal_eps(a, b))
        equal_eps(a, p);
    return orientation_eps(a, b, p) == 0 && greater_eq_eps((b-a)*(p-a), 0.0) && greater_eq_eps((a-b)*(p-b), 0.0);
}

#endif	/* _COMPARE_FLOAT_H */

