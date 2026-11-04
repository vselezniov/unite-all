#ifndef _UNITE_ALL_H
#define	_UNITE_ALL_H

#include "LEDA/incl/LEDA/dictionary.h"


typedef unsigned char BYTE;

// Индкексы в массивах
#define LOWER           0   // Нижняя сторона
#define UPPER           1   // Верхняя сторона

// Флаги ребер
#define LOWER_EDGE      0   // Ребро является нижней стороной трапеции
#define UPPER_EDGE      1   // ... верхней стороной (иначе нижней)
#define RIGHT_TRIANGLE  2   // Ребро принадлежит треугольнику с основанием слева и пиком справа


/*******************************************************************************
 * Реализует основные функции класса |handle|, являющегося всегда не равным NULL
 * указателем на |rep|, ведущим подсчет ссылок
 */
#define IMPLEMENT_HANDLE(handle, rep)        \
    handle() {PTR = new rep();} \
    handle(const handle &h) : handle_base(h) {}     \
    ~handle() {}                                    \
    handle& operator=(const handle &h)              \
    { handle_base::operator=(h); return *this;}     \
    rep* ptr() const; /*{ return (rep*)PTR; }*/   \
    rep* operator->() const { return (rep*)PTR; }


/*********************************************
 *  Трапеция с основаниями параллельными оси Y
 */
struct y_trapezoid
{
    segment edges[2];        // Нижняя и верхняя стороны

    y_trapezoid()
    {}

    y_trapezoid(const segment &lower_edge, const segment &upper_edge)
    { edges[LOWER] = lower_edge; edges[UPPER] = upper_edge; }

    // Приводит трапецию в нормальный вид
    void normalize();
};

/********************************
 * Направленный отрезок с флагами
 */
class flag_segment_rep : public segment_rep
{
    friend class flag_segment;

    int flags;
    
public:
    flag_segment_rep(const point &p, const point &q, int flags)
        : segment_rep(p, q)
    { this->flags = flags; }
    ~flag_segment_rep() {}
};

class flag_segment : public segment
{
    flag_segment_rep* ptr() const { return (flag_segment_rep*)PTR; }

public :

    flag_segment()
        : segment(new flag_segment_rep(point(0,0), point(0,0), 0))
    { }

    flag_segment(const point &p, const point &q, int flags)
        : segment(new flag_segment_rep(p, q, flags))
    { }

    flag_segment(const segment &s, int flags)
        : segment(new flag_segment_rep(s.start(), s.end(), flags))
    { }

    flag_segment(const flag_segment &fs)
        : segment(fs)
    { }

    ~flag_segment() {}

    flag_segment& operator=(const flag_segment& fs)
    { handle_base::operator=(fs); return *this;}

    int flags() const
    { return ptr()->flags; }
};

/********************************************
 * Направленный отрезок с изменяемыми флагами
 */
struct flag_segment2_rep : public handle_rep
{
    segment s;
    int flags;

    flag_segment2_rep()
    { flags = 0; }

    flag_segment2_rep(const segment &s, int flags)
    {
        this->s = s;
        this->flags = flags;
    }

    ~flag_segment2_rep() {};
};

/********************************************************************************
 * Всегда не равный NULL указатель на |flag_segment2_rep|, ведущий подсчет ссылок
 */
class flag_segment2 : public handle_base
{
public:

    IMPLEMENT_HANDLE(flag_segment2, flag_segment2_rep)

    flag_segment2(const segment &s, int flags = 0)
    { PTR = new flag_segment2_rep(s, flags); }

    flag_segment2(const point &start, const point &end, int flags = 0)
    { PTR = new flag_segment2_rep(segment(start, end), flags); }
};

inline int compare(const flag_segment2 &, const flag_segment2 &)
{
    throw "compare(flag_segment2 &, flag_segment2 &): use sweep_flag_segment2_cmp(...) instead";
}

/******************************************************************************
 * Словарь координат, хранящий значения, отличающиеся не менее чем на |EPSILON|
 */
class coord_dict_t
{
    dictionary<double, int> x_dict;
    dictionary<double, int> y_dict;

    // Отыскивает координату в словаре,
    // если такой координаты нет, добавляет её
    double find(const double &v, dictionary<double, int> &dict);

public:

    // Значение, которое больше координат любой точки
    double Infinity;

    coord_dict_t()
    : x_dict(compare_eps), y_dict(compare_eps) {
        Infinity = 1;
    }

    // Отыскивает координату |x| в словаре,
    // если такой координаты нет, добавляет её
    double find_x(const double &x) {
        return find(x, x_dict);
    }

    // Отыскивает координату |y| в словаре,
    // если такой координаты нет, добавляет её
    double find_y(const double &y) {
        return find(y, y_dict);
    }

    // Очищает словарь координат
    void clear()
    {
        x_dict.clear();
        y_dict.clear();
        Infinity = 1;
    }
};


/*
 * Вычисляет объединение многоугольников
 */
gen_polygon unite_all(list<gen_polygon> &poly_list);

/*
 * Формирует набор определенно ориентированных ребер границы многоугольника
 * из набора произвольно ориентированных замкнутых контуров, составляющих границу
 * (это нужно если контуры многоугольника могут не иметь определенной ориентации,
 * т.е. если цепь ребер может сама себя пересекать)
 */
void oriented_boundary(
    /*IN*/      gen_polygon &poly,                  // Входной многоугольник
    /*OUT*/     list<segment> &poly_boundary,       // Выходной набор положительно ориентированных ребер границы
    /*IN, OUT*/ coord_dict_t &coord_dict            // Словарь координат, хранящий значения, отличающиеся не менее чем на |EPSILON|
);

/*
 * Разбивает многоугольник на трапеции с основаниями, параллельными оси Y,
 * и сторонами, принадлежащими границе исходного многоугольника
 */
void trapezoidalize(
    /*IN*/      list<segment> boundary_edges,       // Набор положительно ориентированных отрезков границы многоугольника
    /*OUT*/     list<y_trapezoid> &trapezoids,      // Выходной набор трапеций
    /*IN, OUT*/ coord_dict_t &coord_dict            // Словарь координат, хранящий значения, отличающиеся не менее чем на |EPSILON|
);

/*
 *  Вычисляет невертикальные отрезки границы объединения трапеций
 */
void unite_y_trapezoids(const list<y_trapezoid> &trapezoids,
        list<segment> &uborder_segms, coord_dict_t &coord_dict);

/*
 * Востанавливает вертикальные ребра границы многоугольника по остальным ребрам
 */
void restore_v_edges(list<segment> &boundary_edges, coord_dict_t &coord_dict);

/*
 * Разбивает набор отрезков границы многоугольной области на подмножества
 * по областям связности
 */
void isolate_connected_areas(list<segment> &boundary_edges,
    list< list<segment>* > &connected_areas, coord_dict_t &coord_dict);

/*
 * Формирует замкнутые контуры одного региона связности,
 * добавляет их в набор контуров многоугольной области
 */
void form_contours(list<segment> &connected_edges, list<polygon> &out_contours);

#endif	/* _UNITE_ALL_H */
