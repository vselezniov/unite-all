#ifndef _SWEEP_UNION_H
#define	_SWEEP_UNION_H

#include "LEDA/segment.h"

typedef unsigned char BYTE;

// Индкексы в массивах
#define LOWER           0   // Нижняя сторона
#define UPPER           1   // Верхняя сторона

// Флаги ребер
#define LOWER_EDGE      0   // Ребро является нижней стороной трапеции
#define UPPER_EDGE      1   // ... верхней стороной (иначе нижней)
#define RIGHT_TRIANGLE  2   // Ребро принадлежит треугольнику с основанием слева и пиком справа

// Реализует основные функции класса |handle|, являющегося всегда не равным NULL
// указателем на |rep|, ведущим подсчет ссылок
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
    leda_segment edges[2];        // Нижняя и верхняя стороны

    y_trapezoid()
    {}

    y_trapezoid(const leda_segment &lower_edge, const leda_segment &upper_edge)
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
    flag_segment_rep(const leda_point &p, const leda_point &q, int flags)
        : segment_rep(p, q)
    { this->flags = flags; }
    ~flag_segment_rep() {}
};

class flag_segment : public leda_segment
{
    flag_segment_rep* ptr() const { return (flag_segment_rep*)PTR; }

public :

    flag_segment()
        : leda_segment(new flag_segment_rep(leda_point(0,0), leda_point(0,0), 0))
    { }

    flag_segment(const leda_point &p, const leda_point &q, int flags)
        : leda_segment(new flag_segment_rep(p, q, flags))
    { }

    flag_segment(const leda_segment &s, int flags)
        : leda_segment(new flag_segment_rep(s.start(), s.end(), flags))
    { }

    flag_segment(const flag_segment &fs)
        : leda_segment(fs)
    { }

    ~flag_segment() {}

    flag_segment& operator=(const flag_segment& fs)
    { handle_base::operator=(fs); return *this;}

    int flags() const
    { return ptr()->flags; }
};

struct Polygon;

/*
 * Разбивает многоугольник на трапеции с основаниями, параллельными оси Y,
 * и сторонами, принадлежащими границе исходного многоугольника
 */
void trapezoidalize(const Polygon &poly, leda_list<y_trapezoid> &trapezoids);

/*
 *  Вычисляет невертикальные отрезки границы объединения трапеций
 */
void unite_y_trapezoids(const leda_list<y_trapezoid> &trapezoids, leda_list<leda_segment> &uborder_segms);

/*
 * Востанавливает вертикальные ребра границы многоугольника по остальным ребрам
 */
void restore_v_edges(leda_list<leda_segment> &boundary_edges);

/*
 * Разбивает набор отрезков границы многоугольной области на подмножества
 * по областям связности
 */
void isolate_connected_areas(leda_list<leda_segment> &boundary_edges,
    leda_list< leda_list<leda_segment>* > &connected_areas);

/*
 * Формирует замкнутые контуры одного региона связности,
 * добавляет их в набор контуров многоугольной области
 */
void form_contours(leda_list<leda_segment> &connected_edges, Polygon &outPoly);

#endif	/* _SWEEP_UNION_H */