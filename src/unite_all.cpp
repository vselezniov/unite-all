#include <cassert>
#include <cstdio>
#include <cmath>

#include <iostream>

using namespace std;

#include "LEDA/point.h"
#include "LEDA/segment.h"
#include "LEDA/polygon.h"
#include "LEDA/gen_polygon.h"

int compare(void *const &p1, void *const &p2);

#include "LEDA/list.h"
#include "LEDA/dictionary.h"
#include "LEDA/sortseq.h"
#include "LEDA/p_queue.h"
#include "LEDA/map2.h"

#include "float_compare.h"
#include "unite_all.h"

// Максимальная погрешность при вычислении |double|
const double EPSILON        = 1E-9;
const double EPS_RECIPROC   = 1E+9;		// 1/EPSILON

/***************************************
 * Вычисляет объединение многоугольников
 */
gen_polygon unite_all(list<gen_polygon> &poly_list)
{
    // Словарь координат, хранящий значения, отличающиеся не менее чем на |EPSILON|
    coord_dict_t coord_dict;

    // Разбиваем исходные многоугольники на трапеции
    list<y_trapezoid> trapezoids;
    gen_polygon poly;
    forall(poly, poly_list)
    {
        // Получаем набор положительно ориентированных отреков границы многоугольника
        list<segment> boundary_edges;
        oriented_boundary(poly, boundary_edges, coord_dict);

        // Получаем трапеции многоугольника
        trapezoidalize(boundary_edges, trapezoids, coord_dict);
    }

    // Получаем невертикальные отрезки границы объединения трапеций
    list<segment> union_edges;
    unite_y_trapezoids(trapezoids, union_edges, coord_dict);

    // Восстанавливаем вертикальные ребра границы многоугольника по остальным ребрам
    restore_v_edges(union_edges, coord_dict);

    // Разбиваем набор ребер границы объединения на подмножества по областям
    // связности
    list< list<segment> *> connected_areas;
    isolate_connected_areas(union_edges, connected_areas, coord_dict);

    // Формируем замкнутые контуры из ребер границ областей связности
    list<polygon> out_contours;
    list<segment> *p_connected_edges;
    forall(p_connected_edges, connected_areas)
    {
        form_contours(*p_connected_edges, out_contours);
        delete p_connected_edges;
    }
    connected_areas.clear();

    gen_polygon out_poly = gen_polygon(out_contours, gen_polygon::NO_CHECK);
    return out_poly;
}

/************************************************************************
 * Определяет порядок следования ребер, пересекаемых заметающей прямой
 * при разбиение многоугольника на трапеции
 */
inline int sweep_cmp(const segment &s1, const segment &s2, const point &pnt_event)
{
    /* Предусловие:
     * либо |pnt_event| совпадает с началом одного из отрезков |s1| и |s2|
     * (это означает, что сравнение происходит при вставке нового отрезка
     * в основной список);
     */

    // Накладывающиеся сегменты считаются равными

    if (identical(s1, s2)) return 0;

    int r = 0;

    if (identical(pnt_event, s1.source()))
        // s1 -- новый добавляемый отрезок
        r = half_plane_eps(s2, pnt_event);
    else if (identical(pnt_event, s2.source()))
        // s2 -- новый добавляемый отрезок
        r = -half_plane_eps(s1, pnt_event);
    else
        throw "sweep_cmp(s1,s2): none of the compared edges"
                " starts in the current event-point";

    if (r || s1.is_trivial() || s2.is_trivial())
        return r;

    r = half_plane_eps(s2, s1.target());

    return r;
}

class flag_segment_cmp : public leda_cmp_base<flag_segment>
{
    point pnt_event;  // Текущая точка-событие

public:

    flag_segment_cmp(const point &p) : pnt_event(p)
    { }

    void set_event(const point &p)
    { pnt_event = p; }

    int operator()(const flag_segment &s1, const flag_segment &s2) const
    {
        return sweep_cmp(s1, s2, pnt_event);
    }
};

class flag_segment2_cmp : public leda_cmp_base<flag_segment2>
{
    point pnt_event;  // Текущая точка-событие

public:

    flag_segment2_cmp(const point &p) : pnt_event(p)
    { }

    void set_event(const point &p)
    { pnt_event = p; }

    int operator()(const flag_segment2 &fs1, const flag_segment2 &fs2) const
    {
        return sweep_cmp(fs1->s, fs2->s, pnt_event);
    }
};


/*
 * Вычисляет точку возможного пересечения между заданным ребром и ребром,
 * следующим за ним в наборе |crossed_edges|
 */
inline void compute_intersection(
    /*IN*/      seq_item sit0,                                          // Указатель внутрь |crossed_edges|,
    /*IN*/      sortseq<flag_segment2, seq_item> &crossed_edges,   // Упорядоченная последовательность ребер, пересекаемая заметающей прямой
    /*IN, OUT*/ sortseq<point, seq_item> &event_points,       // Набор точек-событий
    /*IN, OUT*/ coord_dict_t &coord_dict,                               // Словарь координат, хранящий значения, отличающиеся не менее чем на |EPSILON|
    /*IN*/      const point &pnt_event                             // Текущая точка-событие
    )
{
    seq_item sit1 = crossed_edges.succ(sit0);
    segment s0 = crossed_edges.key(sit0)->s;
    segment s1 = crossed_edges.key(sit1)->s;

    int hp1, hp0;

    if ((hp1 = half_plane_eps(s0, s1.target())) <= 0
            && (hp0 = half_plane_eps(s1, s0.target())) >= 0)
    {
        if (hp1 == 0)
            // Конец ребра |s1| лежит на ребре |s0|
            crossed_edges.change_inf(sit0, event_points.lookup(s1.target()));
        else if (hp0 == 0)
            // Конец ребра |s0| лежит на ребре |s1|
            crossed_edges.change_inf(sit0, event_points.lookup(s0.target()));
        else
        {
            // Ребра пересекаются во внутренней точке
            point q;
            s0.intersection_of_lines(s1, q);
            q = point(coord_dict.find_x(q.xcoord()), coord_dict.find_y(q.ycoord()));

            // Проверки связаны с тем, что у почти параллельных близкорасположенных
            // друг к другу ребер вычисленная точка пересеченая может
            // выйти за дипазон
            if(  !( compare(pnt_event, q) < 0 && compare(q, s0.end()) < 0 && compare(q, s1.end()) < 0 )  )
                throw "oriented_boundary(...) compute_intersection(...): computed intersection point is out of range";

            crossed_edges.change_inf(sit0, event_points.insert(q, sit0));
        }
    }
}

/*************************************************************************
 * Формирует набор определенно ориентированных ребер границы многоугольника
 * из набора произвольно ориентированных замкнутых контуров, составляющих границу
 * (это нужно если контуры многоугольника могут не иметь определенной ориентации,
 * т.е. если цепь ребер может сама себя пересекать)
 */
void oriented_boundary(
    /*IN*/      gen_polygon &poly,                  // Входной многоугольник
    /*OUT*/     list<segment> &poly_boundary,       // Выходной набор положительно ориентированных ребер границы
    /*IN, OUT*/ coord_dict_t &coord_dict            // Словарь координат, хранящий значения, отличающиеся не менее чем на |EPSILON|
    )
{
    point pnt_event;                   // Очередная точка-событие
    flag_segment2_cmp cmp(pnt_event);       // Определяет порядок следования пересекаемых ребер

    // Набор точек-событий
    // Для каждой точки дополнительно хранится одно из ребер, оканчивающиеся в ней
    sortseq<point, seq_item> event_points;

    // Очередь ребер на добавление
    p_queue<point, segment> edge_queue;

    // Набор ребер, пересекаемых заметающей прямой,
    // Для каждого ребра дополнительно хранится указатель на точку-событие
    // пересечения со  следущим по порядку ребром в последовательности
    sortseq<flag_segment2, seq_item> crossed_edges(cmp);

    // Таблица, сопоставляющая паре пересекающихся ребер
    // точку-событие, в которой они пересекаются
    map2<flag_segment2, flag_segment2, seq_item> inter_dic;

    // Использовать ли таблицу |inter_dic|
    bool use_optimization = false;

    // Контуры многоугольника
    const list<polygon> &contours = poly.polygons();

    ///////////////////////////////////////////////////////////
    // Инициализируем |event_points|, |edge_queue| и |coord_dict.Infinity|

    // Цикл по контурам многоугольника
    polygon cont;
    forall(cont, contours)
    {
        const list<point> &vertices = cont.vertices();
        int n = vertices.size();

        if (!n)
            continue;

        // Первой точке замкнутого контура предшествует последняя
        point vert = vertices.back();
        vert = point(coord_dict.find_x(vert.xcoord()), coord_dict.find_y(vert.ycoord()));

        seq_item prev_it = event_points.insert(vert, seq_item(nil));
        point prev_pnt = event_points.key(prev_it);

        // Цикл по вершинам контура
        forall(vert, vertices)
        {
            vert = point(coord_dict.find_x(vert.xcoord()), coord_dict.find_y(vert.ycoord()));

            double x = vert.xcoord();
            double y = vert.ycoord();

            if (x < 0) x = -x;
            if (y < 0) y = -y;

            while (x >= coord_dict.Infinity || y >= coord_dict.Infinity)
                coord_dict.Infinity *= 2;

            seq_item it = event_points.insert(vert, seq_item(nil));
            point p = prev_pnt, q = event_points.key(it);

            prev_pnt = q;       // далее возможно будет выполнен "leda_swap(p, q)"

            if(prev_it != it)   // Нас интересуют ребра ненулевой длины
            {
                if(0 < compare(p, q))
                    leda_swap(p, q);

                edge_queue.insert(p, segment(p, q));
            }

            prev_it = it;
        }
    }

    double Infinity = coord_dict.Infinity * 2;      // На всякий случай (для сравнения с использованием epsilon)

    // Вставляем в |crossed_edges| нижний и верхний ограничители
    // для упрощения реализации
    // (флаги этих ограничителей могут использоваться)
    flag_segment2 lower_sentinel(segment(-Infinity, -Infinity, Infinity, -Infinity), UPPER_EDGE);
    flag_segment2 upper_sentinel(segment(-Infinity, Infinity, Infinity, Infinity), LOWER_EDGE);

    // Обработка событий начинается с левой нижней точки
    pnt_event = lower_sentinel->s.start();
    cmp.set_event(pnt_event);

    crossed_edges.insert(upper_sentinel, nil);
    crossed_edges.insert(lower_sentinel, nil);

    // Добавляем ограничитель в |edge_queue|
    point pstop(Infinity, Infinity);
    edge_queue.insert(pstop, segment());

    flag_segment2 fs;   // Одно из ребер в |crossed_edges|
    segment s;     // Ребро результат

    // Очередное добавляемое ребро
    segment next_edge = edge_queue.inf(edge_queue.find_min());

    // Набор указателей внутрь |crossed_edges| на добавленные в текущей
    // точке-событии ребра
    // (используется для выявления противоречий)
    dictionary<void *, int> inserted_dict;

    ////////////////
    // Основной цикл
    while (!event_points.empty())
    {
        // Переходим к следущей точке-событию
        seq_item event_it = event_points.min();
        pnt_event = event_points.key(event_it);
        cmp.set_event(pnt_event);

        //////////////////////////////////////////////////////////////////
        // Обрабатываем ребра заканчивающиеся в текущей точке и проходящие
        // через неё

        // Получаем одно из ребер, проходящее через текущую точку-событие
        seq_item sit = event_points.inf(event_it);

        if (sit == nil)
            // В |crossed_edges| нет ребер оканчивающихся в текущей точке
            // Получаем ребро, проходжящее через эту точку
            sit = crossed_edges.lookup(flag_segment2(pnt_event, pnt_event, 0));

        seq_item sit_pred = nil;    // Ребро, предшествующее ребрам, содержащим текущую точку
        seq_item sit_succ = nil;    // Ребро, следующее за ними
        seq_item sit_pred_succ;     // Ребро, следовавшее за |sit_pred| до пересечения

        if (sit != 0)
        {
            // Есть ребра, содержащие текущую точку
            // |sit| -- указывает на одно из них

            // Сначала доходим до |sit_pred|
            while (crossed_edges.inf(sit) == event_it)
                sit = crossed_edges.succ(sit);

            sit_succ = crossed_edges.succ(sit);
            seq_item sit_last = sit;

            // У всех ребер кроме последнего, пересекающихся в текущей точке
            // дополнительное поле указывает на текущую точку-событие
            if (use_optimization)
            {
                seq_item xit = crossed_edges.inf(sit_last);
                if (xit)
                {
                    flag_segment2 fs1 = crossed_edges.key(sit_last);
                    flag_segment2 fs2 = crossed_edges.key(sit_succ);
                    inter_dic(fs1, fs2) = xit;
                }
            }

            // Спускаемся вниз до |sit_pred|
            // Удаляем закончившиеся ребра, разбиваем продолжающиеся,
            // Помещаем ребра в результат,
            // Обнуляем дополнительную информацию у продолжающихся ребер
            do
            {
                fs = crossed_edges.key(sit);

                if(identical(fs->s.end(), pnt_event))
                {
                    // Закончившееся ребро
                    // Запоминаем ребро результата
                    if(fs->flags == LOWER_EDGE)
                        s = fs->s;
                    else
                        s =fs->s.reverse();

                    // Удаляем ребро, продвигаемся дальше
                    seq_item it = crossed_edges.pred(sit);
                    crossed_edges.del_item(sit);
                    sit = it;
                }
                else
                {
                    // Продолжающееся ребро
                    // Запоминаем ребро результата
                    if(fs->flags == LOWER_EDGE)
                        s = segment(fs->s.start(), pnt_event);
                    else
                        s = segment(pnt_event, fs->s.start());

                    // Разбиваем ребро
                    fs->s = segment(pnt_event, fs->s.end());
                    fs->flags = -1;    // Ориентацию продолжающегося участка нужно выяснить
                    inserted_dict.insert(sit, 0);

                    // Обнуляем дополнительную информацию
                    crossed_edges.change_inf(sit, nil);

                    // Продвигаемся дальше
                    sit = crossed_edges.pred(sit);
                }

                // Добавляем ребро в результат
                poly_boundary.append(s);
                
            } while (crossed_edges.inf(sit) == event_it);

            // |sit_last| уже может быть удален

            sit_pred = sit;
            seq_item sit_first = crossed_edges.succ(sit_pred);
            sit_pred_succ = sit_first;

            sit = sit_first;

            // Меняем порядок следования пересекшихся ребер на обратный
            if (sit_first != sit_succ)
                crossed_edges.reverse_items(crossed_edges.succ(sit_pred), crossed_edges.pred(sit_succ));
        }

        /////////////////////////////////////////////////////////////////
        // Вставляем в |crossed_edge| ребра, начинающиеся в текущей точке
        while (identical(pnt_event, next_edge.source()))
        {
            //cout << next_edge << endl;
            
            // Получаем позицию в |crossed_edges| для нового ребра
            seq_item s_sit = crossed_edges.locate(flag_segment2(next_edge));
            seq_item p_sit = crossed_edges.pred(s_sit);

            fs = crossed_edges.key(s_sit);
            if (half_plane_eps(fs->s, next_edge.start()) == 0 &&
                    half_plane_eps(fs->s, next_edge.end()) == 0)
            {
                // Найдена пара налагающихся ребер

                // Этот блок понадобился для "разрезанных" многоугольников VDRC,
                // в которых дырка входит в основной контур и связана с внешними
                // ребрами парой налагающихся пебер

//                if (fs->flags() == next_edge.flags())
//                    throw "trapezoidalize(...): overlapped edges with the same"
//                    " orientation found";

                // Определяем неперекрывшуюся часть одного из ребер
                segment dif_edge;
                int c = compare(fs->s.end(), next_edge.end());
                if (c < 0)
                    // Конец ребра |next_edge| не накрыт ребром |fs|
                    dif_edge = segment(fs->s.end(), next_edge.end());
                else if(c > 0)
                    // Конец ребра |fs| не накрыт ребром |next_edge|
                    dif_edge = segment(next_edge.end(), fs->s.end());
                else
                    // Ребра заканчиваются в одной точке
                    ;

                // Добавляем |dif_edge| в очередь на добавление
                if (!dif_edge.is_trivial())
                    edge_queue.insert(dif_edge.source(), dif_edge);

                // Удаляем ребро |s_sit| из |crossed_edges|, потому что
                // следующая его часть не является границей многоугольника
                crossed_edges.del_item(s_sit);

                // Удаляем ссылку на |s_sit| из |event_points|
                seq_item end_it = event_points.lookup(fs->s.end());
                if (event_points.inf(end_it) == s_sit)
                    event_points.change_inf(end_it, nil);

                // Удаляем |s_sit| из |inserted_dict|
                dic_item dic_it = inserted_dict.lookup(s_sit);
                assert(dic_it != nil);
                inserted_dict.del_item(dic_it);
            }
            else
            {
                // Всё нормально, ребро не наложилось
                flag_segment2 new_fs = flag_segment2(next_edge, -1);
                sit = crossed_edges.insert_at(s_sit, new_fs, nil);

                // Проверяем что невозникло противоречий
                fs = crossed_edges.key(sit);
                if(!identical(new_fs, fs))
                    throw "oriented_boundary(...): inserted edge was misplaced";

                inserted_dict.insert(sit, 0);
                event_points.insert(next_edge.end(), sit);

                // Если до этого не было ребер, содержащих текущую точку
                // инициализируем |sit_succ|, |sit_pred| и |sit_pred_succ|
                if(sit_pred == nil)
                {
                    sit_succ = s_sit;
                    sit_pred = p_sit;
                    sit_pred_succ = sit_succ;
                }
            }

            edge_queue.del_min();
            next_edge = edge_queue.inf(edge_queue.find_min());
        }

        /////////////////////////////////////////////////////////////////
        // Проходимся по ребрам, начинающимся в текущей точке
        // и выясняем их ориентацию

        if(sit_pred != nil)
        {
            // Получаем ориентацию предшествующего ребра
            int pred_orient = crossed_edges.key(sit_pred)->flags;
            assert(0 <= pred_orient);

            sit = sit_pred;
            while( (sit = crossed_edges.succ(sit)) != sit_succ )
            {
                fs = crossed_edges.key(sit);
                //cout << fs->s << endl;
                assert(fs->flags == -1);

                // Определяем ориентацию ребра
                fs->flags = (pred_orient == UPPER_EDGE ? LOWER_EDGE : UPPER_EDGE);
                pred_orient = fs->flags;

                // Удаляем |sit| из |inserted_dict|
                dic_item dic_it = inserted_dict.lookup(sit);
                assert(dic_it != nil);
                inserted_dict.del_item(dic_it);
            }

            // Выполняем проверку на противоречия
            if(!inserted_dict.empty())
                throw "oriented_boundary(...): inserted edge was misplaced";
        }

        ////////////////////////////////////////////////////////////////////
        // Вычисляем возможные пересечения между |sit_pred| и следующим
        // за ним ребром, и между |sit_succ| и предшествующим ему ребром

        if (sit_pred != nil)
        {
            // У ребра |sit_pred| поменялось ребро, следующее за ним
            // в упорядоченном списке |crossed_edges|.
            // Обнуляем указатель на точку пересечения со следующим
            // за |sit_pred| ребром и вычисляем возможные пересечения.

            if (!use_optimization)
            {
                crossed_edges.change_inf(sit_pred, seq_item(nil));
                compute_intersection(sit_pred, crossed_edges, event_points, coord_dict, pnt_event);
                sit = crossed_edges.pred(sit_succ);
                if(sit != sit_pred)
                    compute_intersection(sit, crossed_edges, event_points, coord_dict, pnt_event);
            }
            else
            {
//                    seq_item xit = crossed_edges.inf(sit_pred);
//                    if (xit)
//                    {
//                        flag_segment2 fs1 = crossed_edges.key(sit_pred);
//                        flag_segment2 fs2 = crossed_edges.key(sit_pred_succ); // sit_first
//                        inter_dic(fs1, fs2) = xit;
//                        crossed_edges.change_inf(sit_pred, seq_item(nil));
//                    }
//
//                    compute_intersection_dict(sit_pred);
//                    sit = crossed_edges.pred(sit_succ);
//                    if (sit != sit_pred)
//                        compute_intersection_dict(sit);
            }
        }

        // Удаляем событие
        event_points.del_item(event_it);
    }

    assert(crossed_edges.size() == 2); // Остались только ограничитель
    assert(edge_queue.size() == 1);
}

/*************************************************************************
 * Разбивает многоугольник на трапеции с основаниями, параллельными оси Y,
 * и сторонами, принадлежащими границе исходного многоугольника
 */
void trapezoidalize(
    /*IN*/      list<segment> boundary_edges, // Набор положительно ориентированных отрезков границы многоугольника
    /*OUT*/     list<y_trapezoid> &trapezoids,     // Выходной набор трапеций
    /*IN, OUT*/ coord_dict_t &coord_dict                // Словарь координат, хранящий значения, отличающиеся не менее чем на |EPSILON|
    )
{
    point pnt_event;               // Очередная точка-событие
    flag_segment_cmp cmp(pnt_event);    // Определяет порядок следования пересекаемых ребер

    // Набор точек-событий
    // Для каждой точки дополнительно хранится одно из ребер, оканчивающиеся в ней
    sortseq<point, seq_item> event_points;

    // Очередь ребер на добавление
    p_queue<point, flag_segment> edge_queue;

    // Набор ребер, пересекаемых заметающей прямой,
    // Для каждого ребра дополнительно хранится точка, левее которой
    // трапеции выделены
    sortseq<flag_segment, point> crossed_edges(cmp);

    ///////////////////////////////////////////////////////////
    // Инициализируем |event_points|, |edge_queue|

    {
        segment s;
        forall(s, boundary_edges)
        {
            seq_item it1 = event_points.insert(s.source(), seq_item(nil));
            seq_item it2 = event_points.insert(s.target(), seq_item(nil));

            point p = event_points.key(it1);
            point q = event_points.key(it2);

            // Нас интересуют только не вертикальные ребра
            if (p.xcoord() != q.xcoord())
            {
                int flags = 0;
                if( !(p.xcoord() < q.xcoord()) )
                {
                    leda_swap(p, q);
                    flags = UPPER_EDGE;
                }

                edge_queue.insert(p, flag_segment(p, q, flags));
            }
        }
    }

    double Infinity = coord_dict.Infinity * 2;      // На всякий случай (для сравнения с использованием epsilon)

    // Вставляем в |crossed_edges| нижний и верхний ограничители
    // для упрощения реализации
    // (флаги этих ограничителей могут использоваться)
    flag_segment lower_sentinel(segment(-Infinity, -Infinity, Infinity, -Infinity), UPPER_EDGE);
    flag_segment upper_sentinel(segment(-Infinity,  Infinity, Infinity,  Infinity), LOWER_EDGE);

    // Обработка событий начинается с левой нижней точки
    pnt_event = lower_sentinel.start();
    cmp.set_event(pnt_event);

    crossed_edges.insert(upper_sentinel, upper_sentinel.start());
    crossed_edges.insert(lower_sentinel, lower_sentinel.start());

    // Добавляем ограничитель в |edge_queue|
    point pstop(Infinity, Infinity);
    edge_queue.insert(pstop, flag_segment());

    // Очередное добавляемое ребро
    flag_segment next_edge = edge_queue.inf(edge_queue.find_min());

    // Текущее положение заметающей линии
    double sweep_x = pnt_event.xcoord();

    // Верхняя точка правого основания последней выделенной трапеции
    double trapez_last_y = -Infinity;

    ////////////////
    // Основной цикл
    while (!event_points.empty())
    {
        // Переходим к следущей точке-событию
        seq_item event_it = event_points.min();
        pnt_event = event_points.key(event_it);
        cmp.set_event(pnt_event);

        if (sweep_x != pnt_event.xcoord())
        {
            trapez_last_y = -Infinity;
            sweep_x = pnt_event.xcoord();
        }

        //////////////////////////////////////////////////////////////////
        // Обрабатываем ребра заканчивающиеся в текущей точке и проходящие
        // через неё

        // Получаем одно из ребер, проходящее через текущую точку-событие
        seq_item sit = event_points.inf(event_it);

        if (sit == nil)
            // В |crossed_edges| нет ребер оканчивающихся в текущей точке
            // Получаем либо ребро, проходжящее через эту точку,
            // либо следующее по порядку ребро в |crossed_edges|
            sit = crossed_edges.locate(flag_segment(pnt_event, pnt_event, 0));

        flag_segment s = crossed_edges.key(sit);

        seq_item sit_pred = nil; // Ребро, предшествующее ребрам, содержащим текущую точку
        seq_item sit_succ = nil; // Ребро, следующее за ними

        // Интервал пар ребер, между которыми будем выделять трапеции
        seq_item sit_first = nil;
        seq_item sit_last = nil;
        int n = 0;

        if (identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0)
        {
            // Есть ребра, содержащие текущую точку
            // |sit| -- указывает на одно из них
            n = 1;

            // Находим |sit_pred|
            sit_pred = sit;
            do
            {
                ++n;
                sit_pred = crossed_edges.pred(sit_pred);
                s = crossed_edges.key(sit_pred);
            }
            while (identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0);
            --n;

            // Находим |sit_succ|
            sit_succ = sit;
            do
            {
                ++n;
                sit_succ = crossed_edges.succ(sit_succ);
                s = crossed_edges.key(sit_succ);
            }
            while (identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0);
            --n;

            // Определяем интервал пар ребер, между которыми будем выделять трапеции

            sit_first = crossed_edges.succ(sit_pred);
            sit_last = crossed_edges.pred(sit_succ);

            flag_segment edge_first = crossed_edges.key(sit_first);
            flag_segment edge_last = crossed_edges.key(sit_last);

            // Проверяем находится ли площадь ниже самого нижнего ребра,
            // содержащего текущую точку
            if (edge_first.flags() == UPPER_EDGE)
            {
                if (trapez_last_y < pnt_event.ycoord())
                {
                    // Трапеция ещё не выделена
                    sit_first = sit_pred;
                    ++n;
                }
                else
                {
                    // Трапеция ниже уже была выделена
                    seq_item it = crossed_edges.succ(sit_first);
                    --n;

                    if (identical(pnt_event, edge_first.end()))
                        crossed_edges.del_item(sit_first);

                    sit_first = it;
                }
            }

            // Если |sit_first| == |sit_last|, то |sit_last| уже может быть удален

            // Проверяем находится ли площадь выше самого верхнего ребра,
            // содержащего текущую точку
            if (edge_last.flags() == LOWER_EDGE)
            {
                sit_last = sit_succ;
                ++n;
            }
        }
        else
        {
            // В ещё |crossed_edges| нет ребер, содержащих текущую точку
            // Если это не изолированная точка (т.е. в ней будут добавлены
            // новые ребра) и если левее её площадь, то нужно выделить трапецию
            // между ребром ниже и выше

            if (s.flags() == UPPER_EDGE && trapez_last_y < pnt_event.ycoord()
                    /*&& identical(pnt_event, next_edge.source())*/)
            {
                sit_last = sit_succ = sit;
                sit_pred = sit_first = crossed_edges.pred(sit);
                n = 2;
            }
        }

        ///////////////////////////////////////////////////
        // Выделяем трапеции между |sit_first| и |sit_last|

        if (n % 2 != 0)
            throw "trapezoidalize(...): n % 2 != 0";

        seq_item sit1 = sit_first;
        while (0 < n / 2)
        {
            seq_item sit2 = crossed_edges.succ(sit1);

            flag_segment s1 = crossed_edges.key(sit1);
            flag_segment s2 = crossed_edges.key(sit2);

            if (s1.flags() != LOWER_EDGE || s2.flags() != UPPER_EDGE)
                throw "trapezoidalize(...): discordant orientation of adjacent edges";

            // Получаем точки на ребрах, левее которых трапеции уже выделены
            point left_pnt1 = crossed_edges.inf(sit1);
            point left_pnt2 = crossed_edges.inf(sit2);

            if (left_pnt1.xcoord() != left_pnt2.xcoord())
                throw "trapezoidalize(...): discordant height of assignable trapezoid";

            if (!(left_pnt1.xcoord() < pnt_event.xcoord()))
                throw "trapezoidalize(...): nonpositive height of assignable trapezoid";

            point right_pnt1 = (sit1 == sit_pred ?
                    point(sweep_x, coord_dict.find_y(s1.y_proj(sweep_x))) : pnt_event);
            point right_pnt2 = (sit2 == sit_succ ?
                    point(sweep_x, coord_dict.find_y(s2.y_proj(sweep_x))) : pnt_event);

            // Добавляем трапецию в результат
            y_trapezoid trapez = y_trapezoid(segment(left_pnt1, right_pnt1),
                                             segment(left_pnt2, right_pnt2));
            trapezoids.append(trapez);

            // Запоминаем верхнюю точку правого основания последней выделенной трапеции
            trapez_last_y = right_pnt2.ycoord();

            // Продвигаемся дальше
            seq_item it = crossed_edges.succ(sit2);

            // Удаляем закончившиеся ребра
            // и запоминаем точки, левее которых трапеции веделены, для оставшихся
            if (identical(pnt_event, s1.end()))
                crossed_edges.del_item(sit1);
            else
                crossed_edges.change_inf(sit1, right_pnt1);

            if (identical(pnt_event, s2.end()))
                crossed_edges.del_item(sit2);
            else
                crossed_edges.change_inf(sit2, right_pnt2);

            sit1 = it;
            n -= 2;
        }

        /////////////////////////////////////////////////////////////////
        // Вставляем в |crossed_edge| ребра, начинающиеся в текущей точке
        while (identical(pnt_event, next_edge.source()))
        {
            sit = crossed_edges.insert(next_edge, next_edge.source());
            event_points.insert(next_edge.end(), sit);

            edge_queue.del_min();
            next_edge = edge_queue.inf(edge_queue.find_min());
        }

        // Удаляем событие
        event_points.del_item(event_it);
    }

    assert(crossed_edges.size() == 2); // Остались только ограничитель
    assert(edge_queue.size() == 1);
}


/******************************************************************************/


/*
 * Востанавливает вертикальные ребра границы многоугольника по остальным ребрам
 */
void restore_v_edges(list<segment> &boundary_edges, coord_dict_t &coord_dict)
{
    point pnt_event;                   // Очередная точка-событие
    flag_segment_cmp cmp(pnt_event);  // Определяет порядок следования пересекаемых ребер

    // Набор точек-событий
    // Для каждой точки дополнительно хранится одно из ребер, оканчивающиеся в ней
    sortseq<point, seq_item> event_points;

    // Очередь ребер на добавление
    p_queue<point, flag_segment> edge_queue;

    // Набор ребер, пересекаемых заметающей прямой
    sortseq<flag_segment, int> crossed_edges(cmp);

    ///////////////////////////////////////////////////////////
    // Инициализируем |event_points|, |edge_queue|

    {
        segment s;
        forall(s, boundary_edges)
        {
            seq_item it1 = event_points.insert(s.source(), seq_item(nil));
            seq_item it2 = event_points.insert(s.target(), seq_item(nil));

            if (it1 == it2) continue; // игнорируем ребра нулевой длины

            point p = event_points.key(it1);
            point q = event_points.key(it2);

            int flags = 0;
            if(0 < compare(p, q))
            {
                leda_swap(p, q);
                flags = UPPER_EDGE;
            }

            edge_queue.insert(p, flag_segment(p, q, flags));
        }
    }

    double Infinity = coord_dict.Infinity * 2;      // На всякий случай (для сравнения с использованием epsilon)

    // Вставляем в |crossed_edges| нижний и верхний ограничители
    // для упрощения реализации
    // (флаги этих ограничителей могут использоваться)
    flag_segment lower_sentinel(segment(-Infinity, -Infinity, Infinity, -Infinity), UPPER_EDGE);
    flag_segment upper_sentinel(segment(-Infinity,  Infinity, Infinity,  Infinity), LOWER_EDGE);

    // Обработка событий начинается с левой нижней точки
    pnt_event = lower_sentinel.start();
    cmp.set_event(pnt_event);

    crossed_edges.insert(upper_sentinel, 0);
    crossed_edges.insert(lower_sentinel, 0);

    // Добавляем ограничитель в |edge_queue|
    point pstop(Infinity, Infinity);
    edge_queue.insert( pstop, flag_segment() );

    // Очередное добавляемое ребро
    flag_segment next_seg = edge_queue.inf(edge_queue.find_min());

    // Что находится слева/справа от текущего положения заметающей прямой,
    // до текущей точки-события и после (true -- внутренность, false -- наружность)
    bool left_before, right_before, left_after, right_after;

    // Значения |left_after|, |right_after| предыдущей точки, имеющей то же
    // значение X как у текущей точки
    bool x_prev_left = false, x_prev_right = false;

    // Первое (нижнее) и последнее (верхнее) ребро, содержащее текущую точку
    // и начинающиеся левее
    flag_segment left_edge_first, left_edge_last;
    // Флаг того, что такие ребра есть
    bool left_edge_present;

    // Первое (нижнее) и последнее (верхнее) ребро, содержащее текущую точку
    // и заканчивающееся правее
    flag_segment right_edge_first, right_edge_last;
    // Флаг того, что такие ребра есть
    bool right_edge_present;

    // Предыдущая точка-событие с тем же значением X как у текущей точки
    point x_prev_pnt = pnt_event;

    // Левая точка горизонтального ребра
    point h_edge_left_pnt;

    ////////////////
    // Основной цикл
    while (!event_points.empty())
    {
        // Переходим к следущей точке-событию
        seq_item event_it = event_points.min();
        pnt_event = event_points.key(event_it);
        cmp.set_event(pnt_event);

        if(x_prev_pnt.xcoord() != pnt_event.xcoord())
        {
            if(x_prev_left != x_prev_right)
                throw "restore_v_edges(...): incomplete vertical edge";

            x_prev_pnt = point(pnt_event.xcoord(), -Infinity); // Эта точка находится на нижнем ограничителе(!)
            x_prev_left = x_prev_right = false;
        }

        //////////////////////////////////////////////////////////////////
        // Обрабатываем ребра заканчивающиеся в текущей точке и проходящие
        // через неё

        // Сейчас все ребра в |crossed_edges|, проходящие через текущую точку
        // начинаются левее неё.
        // По этим ребрам узнаем что находится чуть ниже и левее, и выше и левее
        // текущей точки-события

        // Получаем одно из ребер, содержащее текущую точку-событие
        seq_item sit = event_points.inf(event_it);

        if (sit == nil)
            // В |crossed_edges| нет ребер оканчивающихся в текущей точке
            // Получаем ребро, проходжящее через эту точку, либо следующее
            // ребро в |crossed_edges|
            sit = crossed_edges.locate(flag_segment(pnt_event, pnt_event, 0));

        flag_segment s = crossed_edges.key(sit);

        seq_item sit_pred = nil;        // Ребро, предшествующее ребрам, содержащим текущую точку
        seq_item sit_succ = nil;        // Ребро, следующее за ними

        left_edge_present = false;
        right_edge_present = false;

        // Определяем |sit_pred| и |sit_succ|
        // Получаем |left_edge_first| и |last_edge_left|
        // Удаляем закончившиеся ребра
        if(identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0)
        {
            // Есть ребра, содержащие текущую точку и начинающиеся левее
            // |sit| -- указывает на одно из них

            // Спускаемся вниз, до |sit_pred|
            do
            {
                sit = crossed_edges.pred(sit);
                s = crossed_edges.key(sit);
            } while(identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0);

            sit_pred = sit;
            sit = crossed_edges.succ(sit_pred);
            left_edge_first = crossed_edges.key(sit);

            // Поднимаемся вверх до |sit_succ|, удаляя закончившиеся ребра.
            // Запоминаем верхнее ребро
            s = left_edge_first;
            do
            {
                left_edge_last = s;
                if(identical(pnt_event, s.target()))
                {
                    // Закончившиеся ребро
                    seq_item it = crossed_edges.succ(sit);
                    crossed_edges.del_item(sit);
                    sit = it;
                }
                else
                    // Проходящее ребро
                    sit = crossed_edges.succ(sit);
                s = crossed_edges.key(sit);
            } while(identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0);

            sit_succ = sit;
            left_edge_present = true;
        }
        else
        {
            // Нет ребер, содержащих текущую точку и начинающихся левее
            sit_succ = sit;
            sit_pred = crossed_edges.pred(sit);
        }

        //////////////////////////////////////////////////////////////////
        // Вставляем в |crossed_edge| ребра, начинающиеся в текущей точке.
        while(identical(pnt_event, next_seg.source()))
        {
            sit = crossed_edges.insert(next_seg, 0);
            event_points.insert(next_seg.end(), sit);

            edge_queue.del_min();
            next_seg = edge_queue.inf(edge_queue.find_min());
        }

        // Сейчас все ребра в |crossed_edges|, проходящие через текущую точку
        // оканчиваются правее неё.
        // По этим ребрам узнаем что находится чуть ниже и правее, и выше и правее
        // текущей точки-события

        // Получаем |right_edge_first| и |right_edge_last|
        sit = crossed_edges.succ(sit_pred);
        if(sit != sit_succ)
        {
            right_edge_first = crossed_edges.key(sit);
            sit = crossed_edges.pred(sit_succ);
            right_edge_last = crossed_edges.key(sit);
            right_edge_present = true;
        }

        ///////////////////////////////////
        // Отслеживаем горизонтальные ребра

        // Отслеживается ли начавшееся ранее горизонтальное ребро
        bool trace_h_edge = (x_prev_left != x_prev_right);

        // Выясняем |left_before| и |right_before| по предыдущей точке-событию
        // либо по ребру ниже текущей точки (если оно разделяет текущую и
        // предыдущую точки-события)
        // Сравнение ориентации с "0" не надежно, если координаты заданы с
        // погрешностью, поэтому эти значения по возможности будут пересчитаны
        flag_segment edge_pred = crossed_edges.key(sit_pred);
        if( half_plane_eps(edge_pred, x_prev_pnt) >= 0 )
        {
            // Ребро ниже не разделяет текущую и предыдущую точки-события
            left_before = x_prev_left;
            right_before = x_prev_right;
        }
        else
        {
            // Ребро ниже разделяет текущую и предыдущую точки-события
            assert(edge_pred.xcoord1() < pnt_event.xcoord());
            if(trace_h_edge)
                throw "restore_v_edges(...): vertical edge is being traced"
                        " when previous and current event-points are separated"
                        " by another edge";
            left_before = right_before = (edge_pred.flags() == LOWER);
        }

        // Выясняем теперь значения |left_before|, |left_after|,
        // |right_before|, |right_after| по ребрам, содержащим текущую точку
        if(left_edge_present)
        {
            left_before = left_edge_first.flags() == UPPER_EDGE;
            left_after = (left_edge_last.flags() == LOWER);
        }
        else
            left_after = left_before;

        if(right_edge_present)
        {
            right_before = (right_edge_first.flags() == UPPER_EDGE);
            right_after = (right_edge_last.flags() == LOWER);
        }
        else
            right_after = right_before;

        // Убеждаемся что потребность в горизонтальном ребре согласованна
        if( trace_h_edge != (left_before != right_before) )
            throw "restore_v_edges(...): need in vertical edge is discordant";

        if(trace_h_edge)
        {
            // Продолжается горизонтальное ребро
            if( (left_before != x_prev_left) || (right_before != x_prev_right) )
                throw "restore_v_edges(...): discordant orientation of"
                        " continuing vertical edge";

            // Выясняем нужно ли его завершить
            if( (left_before != left_after) || (right_before != right_after) )
            {
                // Завершаем
                point p = h_edge_left_pnt;
                point q = pnt_event;

                if(!left_before)
                    leda_swap(p, q);

                boundary_edges.append(segment(p, q));
                trace_h_edge = false;
            }
        }

        // Выясняем нужно ли начать новое вертикальное ребро
        if( (left_after != right_after) && !trace_h_edge)
        {
            h_edge_left_pnt = pnt_event;
            trace_h_edge = true;
        }

        // Устанавливаем |x_prev_pnt|, |x_prev_left|, |x_prev_right|
        x_prev_pnt = pnt_event;
        x_prev_left = left_after;
        x_prev_right = right_after;

        // Удаляем событие
        event_points.del_item(event_it);
    }

    if(x_prev_left != x_prev_right)
        throw "restore_v_edges(...): incomplete vertical edge";
    
    assert(crossed_edges.size() == 2);  // Остались только ограничитель
    assert(edge_queue.size() == 1);
}

/************************************************************************
 * Структура описывающая ребро и его принадленжность к области связности
 */
struct conn_edge_rep : handle_rep
{
    flag_segment edge;              // Ребро
    list_item parent;               // Указатель на родительское ребро в дереве
                                    // ребер одной области связности
                                    // (корень дерева указывает сам на себя)

    int num;                        // Для корня хранит кол-во ребер в дереве
    list<segment>* plist; // Для корня хранит указатель на список,
                                    // в который нужно поместить все ребра
                                    // в дереве

    conn_edge_rep()
    {
        parent = nil;
        num = 1;
        plist = nil;
    }

    conn_edge_rep(const flag_segment &s)
    {
        edge = s;
        parent = nil;
        num = 1;
        plist = nil;
    }
};

// Всегда не равный NULL указатель на |edge_connectivity_rep|, ведущий подсчет ссылок
class conn_edge_h : public handle_base
{
public:

    IMPLEMENT_HANDLE(conn_edge_h, conn_edge_rep)
    
    conn_edge_h(const flag_segment& s) {PTR = new conn_edge_rep(s);}
};

// Находит корень дерева ребер региона связности, которому принадлежит данное ребро
// Подвешивает все узлы, на пути к корню, непосредственно к корню
list_item find_root(list_item it0, list<conn_edge_h> &conn_edges)
{
    // Поднимаемся вверх по дереву ребер одного региона связности
    list_item parent_it = it0;
    list_item it;
    conn_edge_h conn_edge;
    do
    {
        it = parent_it;
        conn_edge = conn_edges.contents(it);
        parent_it = conn_edge->parent;
    } while(it != parent_it);

    list_item root_it = it;

    // Проходим цепочку от |it0| до |root_it| ещё раз, подвешивая узлы
    // непосредственно к корню
    conn_edge = conn_edges.contents(it0);
    it = conn_edge->parent;
    while(it != root_it)
    {
        conn_edge->parent = root_it;
        conn_edge = conn_edges.contents(it);
        it = conn_edge->parent;
    }

    return root_it;
}

// Объединяет два региона связности в один
void merge(list_item it1, list_item it2, list<conn_edge_h> &conn_edges)
{
    list_item root_it1 = find_root(it1, conn_edges);
    list_item root_it2 = find_root(it2, conn_edges);

    if(root_it1 != root_it2)
    {
        conn_edge_h root1 = conn_edges.contents(root_it1);
        conn_edge_h root2 = conn_edges.contents(root_it2);

        if(root1->num < root2->num)
        {
            leda_swap(root1, root2);
            leda_swap(root_it1, root_it2);
        }

        // В дереве с корнем |root1| не меньше узлов чем в дереве с корнем |root2|
        // Делаем |root2| поддеревом |root1|
        root2->parent = root_it1;
        root1->num += root2->num;
    }
}

/************************************************************************
 * Разбивает набор отрезков границы многоугольной области на подмножества
 * по областям связности
 */
void isolate_connected_areas(
    list<segment> &boundary_edges,
    list< list<segment>* > &connected_areas,
    coord_dict_t &coord_dict
    )
{
    // Список ребер с информацией о принадлежности к региону связности
    list<conn_edge_h> conn_edges;

    point pnt_event;                   // Очередная точка-событие
    flag_segment_cmp cmp(pnt_event);  // Определяет порядок следования пересекаемых ребер

    // Набор точек-событий
    // Для каждой точки дополнительно хранится одно из ребер, оканчивающиеся в ней
    sortseq<point, seq_item> event_points;

    // Очередь ребер на добавление, хранит указатели внутрь |conn_edges|
    p_queue<point, list_item> edge_queue;

    // Набор ребер, пересекаемых заметающей прямой
    // Дополнительно хранит указатель внутрь |conn_edges|
    sortseq<flag_segment, list_item> crossed_edges(cmp);

    ///////////////////////////////////////////////////////////
    // Инициализируем |conn_edges|, |event_points|, |edge_queue|

    {
        segment s;
        forall(s, boundary_edges)
        {
            seq_item it1 = event_points.insert(s.source(), seq_item(nil));
            seq_item it2 = event_points.insert(s.target(), seq_item(nil));

            if (it1 == it2) continue; // игнорируем ребра нулевой длины

            point p = event_points.key(it1);
            point q = event_points.key(it2);

            // Направленные снизу вверх (сверху вниз) вертикальные ребра
            // ведут себя как нижние (верхние) невертикальные ребра

            int flags = 0;
            if(0 < compare(p, q))
            {
                leda_swap(p, q);
                flags = UPPER_EDGE;
            }

            flag_segment fs = flag_segment(p, q, flags);
            conn_edge_h conn_edge = conn_edge_h(fs);

            list_item conn_edge_it = conn_edges.append(conn_edge);
            conn_edge->parent = conn_edge_it;   // Сначала каждое ребро является единственным узлом отдельного дерева
            edge_queue.insert(p, conn_edge_it);
        }
    }

    double Infinity = coord_dict.Infinity * 2;      // На всякий случай (для сравнения с использованием epsilon)

    // Вставляем в |crossed_edges| нижний и верхний ограничители
    // для упрощения реализации
    // (флаги этих ограничителей могут использоваться)
    flag_segment lower_sentinel(segment(-Infinity, -Infinity, Infinity, -Infinity), UPPER_EDGE);
    flag_segment upper_sentinel(segment(-Infinity,  Infinity, Infinity,  Infinity), LOWER_EDGE);

    // Обработка событий начинается с левой нижней точки
    pnt_event = lower_sentinel.start();
    cmp.set_event(pnt_event);

    crossed_edges.insert(upper_sentinel, nil);
    crossed_edges.insert(lower_sentinel, nil);

    // Добавляем ограничитель в |edge_queue|
    point pstop(Infinity, Infinity);
    edge_queue.insert(pstop, nil);

    // Указатели внутрь |conn_edges| на первое (нижнее) и последнее (верхнее)
    // ребро, содержащее текущую точку и начинающиеся левее
    list_item left_first_it, left_last_it;
    // Флаг того, что такие ребра есть
    bool left_present;

    // Указатели внутрь |conn_edges| на первое (нижнее) и последнее (верхнее)
    // ребро, содержащее текущую точку и заканчивающееся правее
    list_item right_first_it, right_last_it;
    // Флаг того, что такие ребра есть
    bool right_present;

    ////////////////
    // Основной цикл
    while (!event_points.empty())
    {
        // Переходим к следущей точке-событию
        seq_item event_it = event_points.min();
        pnt_event = event_points.key(event_it);
        cmp.set_event(pnt_event);

        // Получаем одно из ребер, содержащее текущую точку-событие
        seq_item sit = event_points.inf(event_it);

        if (sit == nil)
            // В |crossed_edges| нет ребер оканчивающихся в текущей точке,
            // получаем ребро, проходжящее через эту точку, или ребро,
            // следующее за ним
            sit = crossed_edges.locate(flag_segment(pnt_event, pnt_event, 0));

        flag_segment s = crossed_edges.key(sit);

        seq_item sit_pred = nil;        // Ребро, предшествующее ребрам, содержащим текущую точку
        seq_item sit_succ = nil;        // Ребро, следующее за ними

        seq_item sit_passing = nil;     // Ребро содержащее текущую точку,
                                        // но начинающееся раньше неё и заканчивающееся позже
        
        // Флаги того, что удаляются ребра, содержащие текущую точку,
        // расположенные ниже/выше |sit_passing|
        bool was_before_passing = false, was_after_passing = false;

        left_first_it = left_last_it = right_first_it = right_last_it = nil;
        left_present = right_present = false;

        if(identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0)
        {
            // Есть ребра, содержащие текущую точку
            // |sit| -- указывает на одно из них

            // Спускаемся вниз, до |sit_pred|
            do
            {
                sit = crossed_edges.pred(sit);
                s = crossed_edges.key(sit);
            } while(identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0);

            sit_pred = sit;
            sit = crossed_edges.succ(sit_pred);
            left_first_it = crossed_edges.inf(sit);

            // Поднимаемся вверх до |sit_succ|
            // Соседние ребра, угол между которыми содержит внутренность
            // многоугольника объединяем в один регион связности
            // Так же объединяем в один регион связности последнее ребро,
            // содержащее текущую точку, и ребро, следующее за ним, если
            // между ними внутренность.
            // Удаляем закочившиеся ребра
            // Запоминаем верхнее ребро
            s = crossed_edges.key(sit);
            do
            {
                left_last_it = crossed_edges.inf(sit);
                
                seq_item sit2 = crossed_edges.succ(sit);
                flag_segment s2 = crossed_edges.key(sit2);

                if(s.flags() == s2.flags())
                    throw "isolate_connected_areas(...): discordant orientation"
                            " of adjacent edges";

                if(s.flags() == LOWER_EDGE)
                    // Объединяем |sit| и |sit2| в один регион связности
                    merge(crossed_edges.inf(sit), crossed_edges.inf(sit2), conn_edges);

                if(identical(pnt_event, s.target()))
                {
                    // Закончившиеся ребро
                    crossed_edges.del_item(sit);
                    if(sit_passing == nil)
                        was_before_passing = true;
                    else
                        was_after_passing = true;
                }
                else
                {
                    // Ребро проходящее дальше
                    if(sit_passing == nil)
                        sit_passing = sit;
                    else
                        throw "isolate_connected_areas(...): crossing or"
                                " overlapping edges found";
                }
                sit = sit2;
                s = s2;
            } while(identical(s.end(), pnt_event) || half_plane_eps(s, pnt_event) == 0);

            sit_succ = sit;
            left_present = true;
        }
        else
        {
            // Пока что нет ребер, содержащих текущую точку
            sit_succ = sit;
            sit_pred = crossed_edges.pred(sit);
        }

        //////////////////////////////////////////////////////////////////
        // Вставляем в |crossed_edge| ребра, начинающиеся в текущей точке.
        pq_item queue_it;
        while( identical(pnt_event, edge_queue.prio(queue_it = edge_queue.find_min())) )
        {
            list_item conn_edge_it = edge_queue.inf(queue_it);
            conn_edge_h conn_edge = conn_edges.contents(conn_edge_it);

            sit = crossed_edges.insert(conn_edge->edge, conn_edge_it);
            event_points.insert(conn_edge->edge.end(), sit);

            edge_queue.del_min();
        }

        //////////////////////////////////////////////////////////////////////
        // Если есть ребро, содержащее текущую точку |pnt_event|, начинающееся
        // раньше и заканчивающееся позже, и если |pnt_event| является для него
        // точкой внутреннего касания, то это ребро нужно разбить,
        // потому что разные участки этого ребра могут пренадлежать границам
        // разных областей всязности

        if(sit_passing != nil)
        {
            s = crossed_edges.key(sit_passing);

            // Если были удалены или добавлены ребра со стороны площади, то
            // |pnt_event| -- точка внутреннего касания для проходящего ребра.

            bool inner_touch = false;

            if(s.flags() == LOWER_EDGE)
                // |sit_passing| -- нижнее ребро
                inner_touch = ( was_after_passing || (crossed_edges.pred(sit_succ) != sit_passing) );
            else
                // |sit_passing| -- верхнее ребро
                inner_touch = ( was_before_passing || (crossed_edges.succ(sit_pred) != sit_passing) );

            if(inner_touch)
            {
                // Разбиваем это ребро
                // Значение |crossed_edges.key(sit_passing)| не меняем
                // Меняем |conn_edges| и |crossed_edges.inf(sit_passing)|

                list_item conn_edge_it = crossed_edges.inf(sit_passing);
                conn_edge_h conn_edge = conn_edges.contents(conn_edge_it);

                flag_segment s1 = flag_segment(conn_edge->edge.source(), pnt_event, conn_edge->edge.flags());
                flag_segment s2 = flag_segment(pnt_event, conn_edge->edge.end(), conn_edge->edge.flags());

                //cout << s1 << "\n" << s2 << endl;

                // В прежний элемент списка |conn_edges| заносим участок ребра до |pnt_event|
                conn_edge->edge = s1;

                // Создаем новый элемент списка |conn_edges| для второго участка ребра
                conn_edge = conn_edge_h(s2);
                conn_edge_it = conn_edges.append(conn_edge);
                conn_edge->parent = conn_edge_it;   // Пока он в отдельном дереве

                // Далее продолжается второй участок ребра
                crossed_edges.change_inf(sit_passing, conn_edge_it);
            }
        }

        ////////////////////////////////////////////////////////////////////////
        // Обрабатываем ребра, содержащие текущую точку и заканчивающиеся правее

        sit = crossed_edges.succ(sit_pred);
        if(sit != sit_succ)
        {
            // Получаем |right_first_it| и |right_last_it|
            right_first_it = crossed_edges.inf(sit);
            sit = crossed_edges.pred(sit_succ);
            right_last_it = crossed_edges.inf(sit);
            right_present = true;

            // Спускаемся до |sit_pred| по ребрам, содержащим текущую точку и
            // продолжающимся далее
            // Соседние ребра, угол между которыми содержит внутренность
            // многоугольника объединяем в один регион связности
            // (можно начинать с предыдущего ребра, потому что оно тоже продолжается
            // за текущее положение заметающей прямой)
            s = crossed_edges.key(sit);
            do
            {
                seq_item sit0 = crossed_edges.pred(sit);
                flag_segment s0 = crossed_edges.key(sit0);

                if(s0.flags() == s.flags())
                    throw "isolate_connected_areas(...): discordant orientation"
                            " of adjacent edges";

                if(s0.flags() == LOWER_EDGE)
                    // Объединяем |sit| и |sit2| в один регион связности
                    merge(crossed_edges.inf(sit0), crossed_edges.inf(sit), conn_edges);

                sit = sit0;
                s = s0;
            } while(sit != sit_pred);
        }

        ///////////////////////////////////////////////////////////////////////////
        // Теперь обрабатываем углы между первыми/последними левыми/правыми ребрами
        if(left_present != right_present)
        {
            // Есть только левые или только правые ребра, содержащие текущую точку
            list_item first_it = (left_present ? left_first_it : right_first_it);
            list_item last_it = (left_present ? left_last_it : right_last_it);

            conn_edge_h first = conn_edges.contents(first_it);
            conn_edge_h last = conn_edges.contents(last_it);

            if(first->edge.flags() == last->edge.flags())
                throw "isolate_connected_areas(...): discordant orientation"
                            " of adjacent edges";

            if(first->edge.flags() == UPPER_EDGE)
                merge(first_it, last_it, conn_edges);
        }
        else if(left_present && right_present)
        {
            // Есть и левые и правые ребра, содержащие текущую точку
            list_item its[2][2] = { {left_first_it, right_first_it},
                                    {left_last_it, right_last_it} };
            for(int i = 0; i < 2; ++i)
            {
                conn_edge_h left = conn_edges.contents(its[i][0]);
                conn_edge_h right = conn_edges.contents(its[i][1]);

                if(left->edge.flags() != right->edge.flags())
                    throw "isolate_connected_areas(...): discordant orientation"
                            " of adjacent edges";
                if( left->edge.flags() == (i == 0 ? UPPER_EDGE : LOWER_EDGE) )
                    merge(its[i][0], its[i][1], conn_edges);
            }
        }
        
        // Удаляем событие
        event_points.del_item(event_it);
    }

    assert(crossed_edges.size() == 2);  // Остались только ограничители
    assert(edge_queue.size() == 1);

    ///////////////////////////////////////////////////////////
    // Теперь все связные регионы разбиты по отдельным деревьям
    // Заносим узлы каждого дерева в отдельный список
    list_item conn_edge_it, root_it;
    conn_edge_h conn_edge, root;
    forall_items(conn_edge_it, conn_edges)
    {
        conn_edge = conn_edges.contents(conn_edge_it);

        // Находим корень
        root_it = find_root(conn_edge_it, conn_edges);
        root = conn_edges.contents(root_it);
        
        // Создаем отдельный список для области связности, если он ещё не создан
        if(root->plist == nil)
        {
           root->plist = new list<segment>;
           connected_areas.append(root->plist);
        }
        
        // Добавляем ребро в список области связности
        point p = conn_edge->edge.source();
        point q = conn_edge->edge.end();
        if(conn_edge->edge.flags() == UPPER_EDGE)
            leda_swap(p, q);

        root->plist->append(segment(p, q));
    }
}

/************************************************************************
 * Определяет порядок следование ребер в списке, отсортированном по
 * начальным точкам ребер
 */
class start_point_cmp : public leda_cmp_base<point>
{
public:

    bool duplicates;    // Считать ли одинаковые точки разными
                        // (чтобы можно было вставлять дубликаты)

    start_point_cmp()
    { duplicates = true; }

    int operator()(const point &p1, const point &p2) const;
};

int start_point_cmp::operator()(const point &p1, const point &p2) const
{
    int r = compare(p1, p2);

    if(duplicates && r == 0)
    {
        // Считаем одинаковые точки разными,
        // чтобы можно было вставлять дубликаты
        //return (ID_Number(s1) - ID_Number(s2));
        return -1;
    }

    return r;
}

/*******************************************************
 * Формирует замкнутые контуры одного региона связности,
 * добавляет их в многоугольную область
 */
void form_contours(list<segment> &connected_edges, list<polygon> &out_contours)
{
    // Сравнивает ребра по начальным точкам
    start_point_cmp cmp;

    // Ребра связной области, отсортированные по начальным точкам
    sortseq<point, segment> sorted_edges(cmp);

    bool main_cont_found = false;            // Найден ли основной контур
    segment s;

    // Формируем |sorted_edges|
    cmp.duplicates = true;
    forall(s, connected_edges)
        sorted_edges.insert(s.start(), s);

    // Создаем место для основного контура
    list_item main_cont_it = out_contours.append(polygon());

    // Цикл по контурам области связности
    cmp.duplicates = false;
    while(!sorted_edges.empty())
    {
        list<point> point_list;     // Точки контура

        // Получаем первое ребро контура
        seq_item it = sorted_edges.first_item();
        s = sorted_edges.inf(it);
        sorted_edges.del_item(it);

        // Запоминаем начальную точку
        point start_pnt = s.start();
        point prev_pnt = start_pnt, curr_pnt = s.end(), next_pnt;

        // Цикл по вершинам контура
        while(curr_pnt != start_pnt)
        {
            // Получаем следущее ребро и следущую точку
            it = sorted_edges.lookup(curr_pnt);
            if(it == nil)
                throw "form_contours(...): nonclosed contour";

            s = sorted_edges.inf(it);
            sorted_edges.del_item(it);
            next_pnt = s.end();

            if( !collinear_eps(prev_pnt, next_pnt, curr_pnt) )
            {
                // Текущая точка не лежит на одной прямой с предыдущей и следующей,
                // заносим её в контур
                point_list.append(curr_pnt);
                prev_pnt = curr_pnt;
                curr_pnt = next_pnt;
            }
            else
                curr_pnt = next_pnt;
        }

        // |start_pnt| ещё не занесена в контур
        if(point_list.size() < 2)
            throw "form_contours(...): degenerate contour";

        if( !collinear(point_list.back(), point_list.front(), start_pnt) )
            point_list.append(start_pnt);

        // Создаем новый контур
        //cout << point_list << endl;
        polygon cont = polygon(point_list, polygon::NO_CHECK);

        // Выясняем является ли контур основным или дыркой
        int cont_orient = cont.orientation();

        if(cont_orient == 1)
        {
            // Основной контур
            // Ставим его на первую позицию
            
            if(main_cont_found)
                throw "form_contours(...): several main contours found in the area"
                        " of connectedness";

            out_contours.assign(main_cont_it, cont);
            main_cont_found = true;
        }
        else
        {
            // Дырка
            // Добвляем её в список контуров
            out_contours.append(cont);
        }
    }
}

// Отыскивает координату в словаре,
// если такой координаты нет, добавляет её
double coord_dict_t::find(const double &v, dictionary<double, int> &dict)
{
    dic_item it = dict.insert(v, 0);
    double ret_v = dict.key(it);
    if (ret_v != v)
    {
        double d = v - ret_v;
        if (d < 0.0)
        {
            d = -d;
            dic_item it0 = dict.pred(it);
            if (it0 != nil)
            {
                double ret_v0 = dict.key(it0);
                double d0 = v - ret_v0;
                assert(0 < d0);
                if (d0 < d)
                    ret_v = ret_v0;
            }
        }
        else
        {
            dic_item it2 = dict.succ(it);
            if (it2 != nil)
            {
                double ret_v2 = dict.key(it2);
                double d2 = ret_v2 - v;
                assert(0 < d2);
                if (d2 < d)
                    ret_v = ret_v2;
            }
        }
    }
    return ret_v;
}

conn_edge_rep* conn_edge_h::ptr() const { return (conn_edge_rep*)PTR; }

flag_segment2_rep* flag_segment2::ptr() const { return (flag_segment2_rep*)PTR; }
