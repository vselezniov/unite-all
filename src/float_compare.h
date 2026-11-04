#ifndef _FLOAT_COMPARE_H
#define	_FLOAT_COMPARE_H

// Максимальная погрешность при вычислении |double|
extern const double EPSILON;
extern const double EPS_RECIPROC;		// 1/EPSILON

// Функции сравнения, учитывающие погрешность
inline int compare_eps(const double &x, const double &y, const double &eps)
{
    if (x + eps <= y)
        return -1;
    else if (x >= y + eps)
        return 1;
    else
        return 0;
}

inline int compare_eps(const double &x, const double &y)
{
	return compare_eps(x, y, EPSILON);
}

// Равно
inline bool equal_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return ( (x - eps < y) && (y < x + eps) );
}

// Больше
inline bool greater_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x >= y + eps;
}

// Меньше или равно
inline bool less_eq_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x < y + eps;
}

// Меньше
inline bool less_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x <= y - eps;
}

// Больше или равно
inline bool greater_eq_eps(const double &x, const double &y, const double &eps = EPSILON)
{
    return x > y - eps;
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

// Определяет в какой полуплоскости лежит точка относительно направленного отрезка

inline int half_plane_eps(const leda_segment& s, const leda_point& p)
{
    double area = s.dy()*(s.xcoord1() - p.xcoord()) - s.dx()*(s.ycoord1() - p.ycoord());
    if (EPS_RECIPROC*area*area < s.sqr_length()*EPSILON)
        return 0;
    else if (area < 0.0)
        return -1;
    else
        return 1;
}

// Лежат ли три точки на одной прямой
inline bool collinear_eps(const leda_point& a, const leda_point& b, const leda_point& c)
{
    double dx = b.xcoord() - a.xcoord();
    double dy = b.ycoord() - a.ycoord();
    double area = dy * (a.xcoord() - c.xcoord()) - dx * (a.ycoord() - c.ycoord());

    return ( EPS_RECIPROC*area*area < (dx*dx + dy*dy)*EPSILON );
}

#endif	/* _FLOAT_COMPARE_H */

