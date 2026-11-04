#include <cassert>
#include <climits>
#include <iostream>

#include "trapez_union.h"
#include "compare_float.h"

#include <LEDA/p_queue.h>
#include <LEDA/map2.h>

using std::cout;
using std::endl;

/************************************************************************
 * Определяет порядок следования ребер, пересекаемых заметающей прямой
 */
class edge_overlap_cmp : public leda_cmp_base<edge_overlap_h>
{
    point pnt_event;  // Текущая точка-событие

public:

    edge_overlap_cmp(const point &p) : pnt_event(p)
    {}

    void set_event(const point &p)
    { pnt_event = p; }

    int operator()(const edge_overlap_h &edge_ovlp1, const edge_overlap_h &edge_ovlp2) const
    {
        segment s1 = edge_ovlp1->s;
        segment s2 = edge_ovlp2->s;

        //cout << s1 << endl;
        //cout << s2 << endl;

        /* Предусловие:
         * либо |pnt_event| совпадает с началом одного из отрезков |s1| и |s2|
         * (это означает, что сравнение происходит при вставке нового отрезка
         * в основной список);
         * либо оба отрезка содержат |pnt_event|,
         * (это означает что сравнение выполняется при подготоке вспомогательного
         * списка добавляемых/удаляемых ребер, проходящих через текущую точку)
         */

        // Накладывающиеся сегменты считаются равными

        if (identical(s1, s2)) return 0;

        int s = 0;

        if (identical(pnt_event, s1.source()))
            // s1 -- новый добавляемый отрезок
            s = orientation_eps(s2, pnt_event);
        else if (identical(pnt_event, s2.source()))
            // s2 -- новый добавляемый отрезок
            s = -orientation_eps(s1, pnt_event);
        else
        {
            // Сравниваются два отрезка, пересекающиеся или оканчивающиеся
            // в текущей точке-событии
            assert(contains_eps(s1, pnt_event) && contains_eps(s2, pnt_event));

            s = orientation_eps(s1, s2.source());
            return s;
        }

        if (s || s1.is_trivial() || s2.is_trivial())
            return s;

        s = orientation_eps(s2, s1.target());

        return s;
    }
};

/*
 * Вычисляет точку возможного пересечения между заданным ребром и ребром,
 * следующим за ним в наборе |crossed_edges|
 */
static void compute_intersection(
    sortseq<point, seq_item>            &event_points,      // Набор точек-событий
    sortseq<edge_overlap_h, seq_item>   &crossed_edges,     // Набор ребер, пересекаемых заметающей прямой
    seq_item edge_it0                                           // Указатель внутрь |crossed_edges|
    )
{
    seq_item edge_it1 = crossed_edges.succ(edge_it0);
    segment s0 = crossed_edges.key(edge_it0)->s;
    segment s1 = crossed_edges.key(edge_it1)->s;

    if (orientation_eps(s0, s1.target()) <= 0 &&
            orientation_eps(s1, s0.target()) >= 0)
    {
        point q;
        s0.intersection_of_lines(s1, q);
        crossed_edges.change_inf(edge_it0, event_points.insert(q, edge_it0));
    }
}

/*
 * Вычисляет точку возможного пересечения между заданным ребром и ребром,
 * следующим за ним в наборе |crossed_edges|
 */
static void compute_intersection(
    sortseq<point, seq_item>            &event_points,      // Набор точек-событий
    sortseq<edge_overlap_h, seq_item>   &crossed_edges,     // Набор ребер, пересекаемых заметающей прямой
    map2<edge_overlap_h, edge_overlap_h, seq_item> &inter_dic,  // Таблица пересечений
    seq_item edge_it0                                           // Указатель внутрь |crossed_edges|
    )
{
    // Given an item |edge_it0| in the Y-structure compute the point of
    // intersection with its successor and (if existing) insert it into
    // the event queue and do all necessary updates.

    seq_item edge_it1 = crossed_edges.succ(edge_it0);

    edge_overlap_h edge_ovlp1 = crossed_edges.key(edge_it0);
    edge_overlap_h edge_ovlp2 = crossed_edges.key(edge_it1);

    segment s0 = edge_ovlp1->s;
    segment s1 = edge_ovlp2->s;

    // |s1| is the successor of |s0| in the Y-structure, hence,
    // |s0| and |s1| intersect right or above of the sweep line
    // iff |(s0.start(),s0.end(),s1.end()| is not a left turn and
    // |(s1.start(),s1.end(),s0.end()| is not a right turn.
    // In this case we intersect the underlying lines

    if (orientation_eps(s0, s1.target()) <= 0 &&
            orientation_eps(s1, s0.target()) >= 0)
    {
        seq_item it = inter_dic.operator()(edge_ovlp1, edge_ovlp2);
        //seq_item it = inter_dic(s0,s1);
        if (it == nil)
        {
            point q;
            s0.intersection_of_lines(s1, q);
            it = event_points.insert(q, edge_it0);
        }
        crossed_edges.change_inf(edge_it0, it);
    }
}

/***********************************************************************
 * Содержит данные и вспомогательные функции, необходимые для вычисления
 * границы объединения трапеций с основаниями, параллельными оси Y
 */
struct y_trapez_sweep_box
{
    // Очередная точка-событие (начало/конец или пересечение ребер)
    point pnt_event;

    // Указатель на очередную точку-событие
    seq_item event_it;
    
    // Определяет порядок следования ребер, пересекаемых заметающей прямой
    edge_overlap_cmp cmp;

    // Набор точек-событий
    // (для каждой точки дополнительно хранится указатель на одно из
    // пересекаемых ребер, проходящих через эту точку)
    sortseq<point, seq_item> event_points;

    // Очередь ребер на добавление
    p_queue<point, flag_segment> edge_queue;

    // Очередь ребер на удаление
    p_queue<point, flag_seq_item> edge_del_queue;

    // Набор ребер трапеций, пересекаемых (вертикальной) заметающей прямой.
    // Налагающиеся ребра представлены одним элементом
    // (дополнительно хранится указатель на точку-событие пересечения со
    // следущим по порядку ребром в последовательности)
    sortseq<edge_overlap_h, seq_item> crossed_edges;

    // Значение, которое больше координат всех точек
    double Infinity;

    // Указывает на элемент не принадлежащий списку |crossed_edges|
    edge_overlap_h tmp_edge_ovlp;

    // Ребро, следующее за пересекающимися в текущей точке
    seq_item edge_it_succ;
    // Ребро, предшествующее пересекающимся
    seq_item edge_it_pred;
    // Ребро, следовавшее за |edge_it_pred| до пересечение
    seq_item edge_it_pred_succ;
    // Одно из ребер, содержащих текущую точку-событие
    seq_item edge_it;

    // Список ребер, удаляемых в текущей-точке событии
    // (ребра правых треугольников в него не заносятся,
    //  они удаляются сразу по извлечению из очереди на удаление)
    list<flag_seq_item> event_del_edges;

    // Подпоследовательность последовательности |crossed_edges|.
    // Каждый элемент предшествует сплошному интервалу ребер,
    // на котором либо все счетчики |inner_cnt| изменились,
    // либо все счетчики |inner_cnt| не изменились
    list<seq_item> inner_intervals;

    // Указатель на список невертикальных ребер объединения входных трапеций
    list<segment> *p_union_edges;

    // Таблица, сопоставляющая паре пересекающихся ребер
    // точку-событие, в которой они пересекаются
    map2<edge_overlap_h, edge_overlap_h, seq_item> inter_dic;

    // Использовать ли таблицу |inter_dic|
    bool use_optimization;

    y_trapez_sweep_box()
    :   cmp(pnt_event),
        event_points(compare_eps),
        edge_queue(compare_eps),
        edge_del_queue(compare_eps),
        crossed_edges(cmp),
        inter_dic(nil)
    {
        use_optimization = false;
        Infinity = 1;
        edge_it_succ = edge_it_pred = edge_it_pred_succ = edge_it = nil;
        p_union_edges = nil;
    }

    // Инициализирует структуры данных
    void init(const list<y_trapezoid> &trapezoids);

    // Устанавливает текущую точку-событие
    void set_event(const point &pnt)
    {
        pnt_event = pnt;
        cmp.set_event(pnt);
    }

    // Вычисляет невертикальные отрезки границы объединения трапеций
    void unite_y_trapezoids(const list<y_trapezoid> &trapezoids, list<segment> &uborder_segms);

    // Формирует список ребер удаляемых в текущей точке
    void prepare_del_edges();

    // Пересчитывает счетчики |inner_cnt| у ребер, пересекающихся в текущей
    // точке-событии, меняет порядок следования пересекшихся ребер на обратный
    void intersect_edges();

    // Вычитатает из счетчиков наложений удаляемые ребра
    void delete_edges();

    // Обновляет флаги произведенных изменений в группе наложившихся ребер
    void update_change_flags(
        seq_item edge_ovlp_it,  // Группа наложившихся ребер
        int egde_orient,        // 0 -- LOWER (edge), 1 -- UPPER
        int change              // (-1) -- удалено, 1 -- добавлено
    );

    // Вставляет ребра, начинающиеся в текущей точке, в список пересекаемых заметающей прямой
    void insert_new_edges();

    // На основании изменений, описанных в |event_changes|
    // Формирует интервалы ребер в последовательности |crossed_edges|,
    // на которых либов все счетчики |inner_cnt| изменились,
    // либо все не изменились.
    // Удаляет из списка |crossed_edges| пустые группы наложившихся ребер
    void append_inner_intervals();

    // Вычисляет возможные пересечения между |edge_it_pred| и следующим за ним ребром,
    // и между |edge_it_succ| и предшествующим ему ребром
    void compute_intersections();

    // Выполняет согласованное обновление счетчиков |inner_cnt|
    // в связи с удалениями/добавлениями ребер в точках с одинаковым X-ом
    void update_inner_counters(double sweep_x);

    // Выставляет значение флагов "граница объединения"
    // и если обнаружен участок границы объединения, добавляет его в |p_union_edges|
    void set_uborder_flags(edge_overlap_h &edge_ovlp, const point &pnt_event);
    void set_uborder_flags(edge_overlap_h &edge_ovlp, double sweep_x);
};

/* Инициализирует структуры данных
 *  Аргументы:
 *      trapezoids  --  список входных трапеций
 */
void y_trapez_sweep_box::init(const list<y_trapezoid> &trapezoids)
{
    Infinity = 1;

    y_trapezoid trapez;
    bool is_right_triangle;     // Является ли трапеция треугольником
                                // с основанием слева и пиком справа

    forall(trapez, trapezoids)
    {
        // Приводим трапецию в нормальный вид
        trapez.normalize();

        is_right_triangle = equal_eps( trapez.edges[LOWER].target(), trapez.edges[UPPER].target() );

        for(int i = 0; i < 2; ++i)
        {
            segment s = trapez.edges[i];

            double x1 = s.xcoord1();
            double y1 = s.ycoord1();
            double x2 = s.xcoord2();
            double y2 = s.ycoord2();
            if (x1 < 0) x1 = -x1;
            if (y1 < 0) y1 = -y1;
            if (x2 < 0) x2 = -x2;
            if (y2 < 0) y2 = -y2;

            while (x1 >= Infinity || y1 >= Infinity ||
                    x2 >= Infinity || y2 >= Infinity)
                Infinity *= 2;

            seq_item it1 = event_points.insert(s.source(), seq_item(nil));
            seq_item it2 = event_points.insert(s.target(), seq_item(nil));

            if (it1 == it2)
                continue; // игнорируем сегменты нулевой длины

            point p = event_points.key(it1);
            point q = event_points.key(it2);
            assert(compare_eps(p, q) < 0);

            int flags = (i == 0 ? LOWER_EDGE : UPPER_EDGE)
                    | (is_right_triangle ? RIGHT_TRIANGLE : 0);
            flag_segment fs1 = flag_segment(p, q, flags);

            edge_queue.insert(fs1.source(), fs1);
        }
    }

    // Вставляем в |crossed_edges| нижний и верхний ограничители
    // для упрощения реализации
    edge_overlap_h lower_sentinel(segment(-Infinity, -Infinity, Infinity, -Infinity));
    edge_overlap_h upper_sentinel(segment(-Infinity,  Infinity, Infinity, Infinity));

    // Ни трапеции ни сегменты их не накрывают (эти значения могут использоваться)
    lower_sentinel->inner_cnt = upper_sentinel->inner_cnt = 0;

    // Обработка событий начинается с левой нижней точки
    set_event(lower_sentinel->s.source());

    crossed_edges.insert(upper_sentinel, seq_item(nil));
    crossed_edges.insert(lower_sentinel, seq_item(nil));

    // Добавляем ограничитель в |edge_queue| и |edge_del_queue|
    point pstop(Infinity, Infinity);
    edge_queue.insert( pstop, flag_segment() );
    edge_del_queue.insert( pstop, flag_seq_item() );
}


/* Формирует список ребер удаляемых в текущей точке */
void y_trapez_sweep_box::prepare_del_edges()
{
    pq_item queue_it;
    while( identical(pnt_event, edge_del_queue.prio(queue_it = edge_del_queue.find_min())) )
    {
        flag_seq_item del_edge = edge_del_queue.inf(queue_it);
        edge_overlap_h edge_ovlp = crossed_edges.key(del_edge.it);
        //cout << edge_ovlp->s;

        bool is_right_triangle = del_edge.flags & RIGHT_TRIANGLE;

        // Обновляем список изменений
        if(!is_right_triangle)
            update_change_flags(del_edge.it, del_edge.flags & UPPER_EDGE, -1);
        else
            // Для правых треугольников ребра, при согласованном изменении
            // счетчика накрытий внутренностями, нижние ребра
            // засчитываются как верхние, а верхние как нижние
            update_change_flags(del_edge.it, 1 - (del_edge.flags & UPPER_EDGE), -1);

        if(is_right_triangle)
            // Ребра правых треугольников удаляем до обработки пересечений
            --(edge_ovlp->edge_cnts[del_edge.flags & UPPER_EDGE]);
        else
            event_del_edges.append(del_edge);

        edge_del_queue.del_min();
    }
}

// Пересчитывает счетчики |inner_cnt| у ребер, пересекающихся в текущей
// точке-событии, меняет порядок пересекшихся ребер на обратный
void y_trapez_sweep_box::intersect_edges()
{
    // edge_it -- указатель на одно из ребер, проходящих через текущую точку-событие
    assert(edge_it != nil);

    // Сначала доходим до |edge_it_succ|

    while (crossed_edges.inf(edge_it) == event_it)
        edge_it = crossed_edges.succ(edge_it);

    edge_it_succ = crossed_edges.succ(edge_it);
    seq_item edge_it_last = edge_it;

    // У всех ребер кроме последнего, пересекающихся в текущей точке
    // дополнительное поле указывает на текущую точку-событие

    if (use_optimization)
    {
        seq_item xit = crossed_edges.inf(edge_it_last);
        if (xit)
        {
            edge_overlap_h edge_ovlp1 = crossed_edges.key(edge_it_last);
            edge_overlap_h edge_ovlp2 = crossed_edges.key(edge_it_succ);
            inter_dic(edge_ovlp1, edge_ovlp2) = xit;
        }
    }

    // У последнего оно уже не актуально, так как порядок сегментов
    // поменяется после точки пересечения
    crossed_edges.change_inf(edge_it_last, seq_item(nil));

    // Спускаемся вниз до |edge_it_pred|, обновляя счетчики |inner_cnt|
    // в связи с пересечением ребер.
    // Обнуляем дополнительную информацию у пересекающихся ребер
    while(crossed_edges.inf(edge_it = crossed_edges.pred(edge_it)) == event_it)
    {
        edge_overlap_h edge_ovlp = crossed_edges.key(edge_it);
        seq_item edge_it2 = edge_it_last;
        do
        {
            edge_overlap_h edge_ovlp2 = crossed_edges.key(edge_it2);

            // |edge_it| ниже |edge_it2| до точки пересечения
            // поэтому |edge_it| входит во внутренности нижних ребер |edge_it2|
            // и выходит из внутренностей верхних ребер.
            // Для |edge_it2| ситуация обратная
            edge_ovlp->inner_cnt += edge_ovlp2->edge_cnts[LOWER] - edge_ovlp2->edge_cnts[UPPER];
            edge_ovlp2->inner_cnt -= edge_ovlp->edge_cnts[LOWER] - edge_ovlp->edge_cnts[UPPER];

            edge_it2 = crossed_edges.pred(edge_it2);
        } while(edge_it2 != edge_it);

        crossed_edges.change_inf(edge_it, seq_item(nil));
    }

    edge_it_pred = edge_it;
    seq_item edge_it_first = crossed_edges.succ(edge_it_pred);
    edge_it_pred_succ = edge_it_first;

    edge_it = edge_it_first;

    // Меняем порядок следования пересекшихся ребер на обратный
    if (edge_it_first != edge_it_last)
        crossed_edges.reverse_items(crossed_edges.succ(edge_it_pred), crossed_edges.pred(edge_it_succ));
}

// Выставляет значение флагов "граница объединения"
// и если обнаружен участок границы объединения, добавляет его в |p_union_edges|
void y_trapez_sweep_box::set_uborder_flags(edge_overlap_h &edge_ovlp, const point &pnt_event)
{
    // Если участок был границей, а стал не границей, то нужно
    // добавить в результат новый найденный отрезок границы и сбросить
    // флаг "граница объединения"

    // Если участок был правой (левой) границей, а стал левой (правой)
    // границей, то нужно добавить в результат новый найденный отрезок
    // границы и выставить новое значение флага

    // Если участок не был границей, а стал границей, то нужно
    // выставить новое значение флага

    bool is_border = (edge_ovlp->inner_cnt == 0 &&
        (0 == edge_ovlp->edge_cnts[LOWER] || 0 == edge_ovlp->edge_cnts[UPPER]) &&
        (0 < edge_ovlp->edge_cnts[LOWER] || 0 < edge_ovlp->edge_cnts[UPPER]) );

    bool was_border = (edge_ovlp->flags & UNION_BORDER);

    bool border_orient_changed = ((0 < edge_ovlp->edge_cnts[UPPER]) != (edge_ovlp->flags & UPPER_BORDER));

    if(  is_border != was_border || (border_orient_changed = (is_border && border_orient_changed) )  )
    {
        if(was_border)
        {
            // Участок был границей, а теперь либо не является границей, либо ориентация
            // границы поменялась
            if( less_eps(edge_ovlp->uborder_since.xcoord(), pnt_event.xcoord()) )
            {
                if(edge_ovlp->flags & UPPER_BORDER)
                    p_union_edges->append(segment(pnt_event, edge_ovlp->uborder_since));
                else
                    p_union_edges->append(segment(edge_ovlp->uborder_since, pnt_event));
            }
        }

        edge_ovlp->flags &= ~BORDER_MASK;

        if(is_border)
        {
            edge_ovlp->flags |= ( UNION_BORDER | (0 < edge_ovlp->edge_cnts[UPPER] ? UPPER_BORDER : LOWER_BORDER) );
            edge_ovlp->uborder_since = pnt_event;
        }
    }
}

// Выставляет значение флагов "граница объединения"
// и если обнаружен участок границы объединения, добавляет его в |p_union_edges|
void y_trapez_sweep_box::set_uborder_flags(edge_overlap_h &edge_ovlp, double sweep_x)
{
    // Если участок был границей, а стал не границей, то нужно
    // добавить в результат новый найденный отрезок границы и сбросить
    // флаг "граница объединения"

    // Если участок был правой (левой) границей, а стал левой (правой)
    // границей, то нужно добавить в результат новый найденный отрезок
    // границы и выставить новое значение флага

    // Если участок не был границей, а стал границей, то нужно
    // выставить новое значение флага

    bool is_border = (edge_ovlp->inner_cnt == 0 &&
        (0 == edge_ovlp->edge_cnts[LOWER] || 0 == edge_ovlp->edge_cnts[UPPER]) &&
        (0 < edge_ovlp->edge_cnts[LOWER] || 0 < edge_ovlp->edge_cnts[UPPER]) );

    bool was_border = (edge_ovlp->flags & UNION_BORDER);

    bool border_orient_changed = ((0 < edge_ovlp->edge_cnts[UPPER]) != (edge_ovlp->flags & UPPER_BORDER));

    if(  is_border != was_border || (border_orient_changed = (is_border && border_orient_changed) )  )
    {
        point pnt = point(sweep_x, edge_ovlp->s.y_proj(sweep_x));

        if(was_border)
        {
            // Участок был границей, а теперь либо не является границей, либо ориентация
            // границы поменялась
            if( less_eps(edge_ovlp->uborder_since.xcoord(), pnt.xcoord()) )
            {
                if(edge_ovlp->flags & UPPER_BORDER)
                    p_union_edges->append(segment(pnt, edge_ovlp->uborder_since));
                else
                    p_union_edges->append(segment(edge_ovlp->uborder_since, pnt));
            }
        }

        edge_ovlp->flags &= ~BORDER_MASK;

        if(is_border)
        {
            edge_ovlp->flags |= ( UNION_BORDER | (0 < edge_ovlp->edge_cnts[UPPER] ? UPPER_BORDER : LOWER_BORDER) );
            edge_ovlp->uborder_since = pnt;
        }
    }
}

// Вычитатает из счетчиков наложений удаляемые ребра
void y_trapez_sweep_box::delete_edges()
{
    flag_seq_item del_edge;
    forall(del_edge, event_del_edges)
    {
        edge_overlap_h edge_ovlp = crossed_edges.key(del_edge.it);

        --(edge_ovlp->edge_cnts[del_edge.flags & UPPER_EDGE]);

        // Обновляем список изменений
        //update_change_list(del_edge.it, del_edge.flags & UPPER_EDGE, -1); // уже сдалано в prepare_del_edges
    }
    event_del_edges.clear();
}

// Обновляет флаги произведенных изменений в группе наложившихся ребер
void y_trapez_sweep_box::update_change_flags(
        seq_item edge_ovlp_it,  // Группа наложившихся ребер
        int egde_orient,        // 0 -- LOWER (edge), 1 -- UPPER
        int change              // (-1) -- удалено, 1 -- добавлено
    )
{
    assert(LOWER <= egde_orient && egde_orient <= UPPER);
    assert(change == -1 || change == 1);
            
    int change_flag = 0;
    if(egde_orient == LOWER)
        change_flag = (change == -1 ? LOWER_DELETED : LOWER_ADDED);
    else
        change_flag = (change == -1 ? UPPER_DELETED : UPPER_ADDED);

    edge_overlap_h edge_ovlp = crossed_edges.key(edge_ovlp_it);
    edge_ovlp->flags |= change_flag;
}

// Вставляет ребра, начинающиеся в текущей точке, в список пересекаемых заметающей прямой
void y_trapez_sweep_box::insert_new_edges()
{
    // Получаем следующее ребро
    pq_item queue_it;
    while( identical(pnt_event, edge_queue.prio(queue_it = edge_queue.find_min())) )
    {
        flag_segment next_edge = edge_queue.inf(queue_it);

        // insert |next_seg| into the Y-structure and associate the
        // corresponding item with the right endpoint of |next_seg| in
        // the X-structure (already present)

        // Находим позицию для нового ребра в последовательности пересекаемых ребер
        tmp_edge_ovlp->s = next_edge;
        seq_item s_edge_it = crossed_edges.locate(tmp_edge_ovlp);
        seq_item p_edge_it = crossed_edges.pred(s_edge_it);

        // Добавляем ребро в список пересекаемых
        edge_overlap_h edge_ovlp = crossed_edges.key(s_edge_it);
        if (orientation_eps(edge_ovlp->s, next_edge.start()) == 0 &&
                orientation_eps(edge_ovlp->s, next_edge.end()) == 0)
        {
            // Ребро наложилось на другое ребро
            if(edge_ovlp->s.xcoord2() < next_edge.xcoord2())
               edge_ovlp->s = segment(edge_ovlp->s.start(), next_edge.end()) ;
            edge_it = s_edge_it;
            s_edge_it = crossed_edges.succ(s_edge_it);
        }
        else
        {
            // Новый пересекаемый прямолинейный участок объединения границ
            edge_it = crossed_edges.insert_at(s_edge_it, tmp_edge_ovlp, seq_item(nil));
            edge_ovlp = tmp_edge_ovlp;
            tmp_edge_ovlp = edge_overlap_h();
        }

        ++edge_ovlp->edge_cnts[next_edge.flags() & UPPER_EDGE];

        // Обновляем список изменений
        update_change_flags(edge_it, next_edge.flags() & UPPER_EDGE, 1);

        // Добавляем ребро в очередь удаления
        edge_del_queue.insert(next_edge.end(), flag_seq_item(edge_it, next_edge.flags()));
        event_points.insert(next_edge.end(), edge_it);

        // Если до этого не было пересекаемых ребер, проходящих через текущую точку
        // инициализируем |edge_it_succ|, |edge_it_pred| и |edge_it_pred_succ|
        if (edge_it_succ == nil)
        {
            edge_it_succ = s_edge_it;
            edge_it_pred = p_edge_it;
            edge_it_pred_succ = edge_it_succ;
        }

        edge_queue.del_min();
    }
}

// На основании флагов произведенных изменений в |edge_overlap_rep::flags|
// Формирует интервалы ребер в последовательности |crossed_edges|,
// на которых либо все счетчики |inner_cnt| изменились,
// либо все не изменились.
// Удаляет из списка |crossed_edges| пустые группы наложившихся ребер
void y_trapez_sweep_box::append_inner_intervals()
{
    if(edge_it_pred == nil)
        return;

    // Получаем последний интервал, уже присутствующий в |inner_intervals|
    seq_item last_interval = nil;
    if( !inner_intervals.empty() )
        last_interval = inner_intervals.back();

    // Цикл по всем ребрам, содержащим текущую точку-событие

    seq_item edge_it = crossed_edges.succ(edge_it_pred);
    while(edge_it != edge_it_succ)
    {
        edge_overlap_h edge_ovlp = crossed_edges.key(edge_it);

        // Устанавливаем флаги "граница объединения"
        // Добавляем ребро  |p_union_edges| при обнаружении участка границы объединения
        set_uborder_flags(edge_ovlp, pnt_event);

        // Получаем флаги, описывающие произведенные изменения
        int change_flags = edge_ovlp->flags & CHANGES_MASK;

        if(change_flags != 0)
        {
            // Выясняем, была ли группа наложившихся ребер добавлена в текущей
            // точке-событии
            assert(edge_ovlp->inner_cnt == INT_MIN || 0 <= edge_ovlp->inner_cnt);
            bool new_edge_ovlp = (edge_ovlp->inner_cnt == INT_MIN);

            // Выясняем, пуста ли группа наложившихся ребер после всех изменений
            assert(0 <= edge_ovlp->edge_cnts[LOWER] && 0 <= edge_ovlp->edge_cnts[UPPER]);
            bool empty_edge_ovlp = (edge_ovlp->edge_cnts[LOWER] == 0 && edge_ovlp->edge_cnts[UPPER] == 0);

            assert(  !empty_edge_ovlp ||
                      ( !(change_flags & (LOWER_ADDED | UPPER_ADDED)) &&
                         (change_flags & (LOWER_DELETED | UPPER_DELETED)) )  );

            // Каждая изменившаяся группа добавляет один или два интеравала
            // обновления счетчиков

            /* Нам нужен интервал, указуемый элементом, предшествующим
             * |edge_it| если:
             *  1) элемент |edge_it| был добавлен в |crossed_edges| в текущей
             *  точке-событии
             *  2) Группа пустая (она удаляется)
             *  2) добавлено верхнее ребро
             *  3) удалено верхнее ребро
             */

            if( new_edge_ovlp || empty_edge_ovlp ||
                    (change_flags & (UPPER_ADDED | UPPER_DELETED)) )
            {
                seq_item prev_edge_it = crossed_edges.pred(edge_it);
                edge_overlap_h prev_edge_ovlp = crossed_edges.key(prev_edge_it);
                //cout << prev_edge_ovlp->s;
                if(last_interval != prev_edge_it)
                {
                    if(prev_edge_ovlp->inner_interval_it != nil)
                        throw "append_inner_intervals(): attemp to append"
                                "one edge twice";

                    prev_edge_ovlp->inner_interval_it =
                            inner_intervals.append(prev_edge_it);
                    last_interval = prev_edge_it;
                }
            }

            /* Нам нужен интервал, указуемый |edge_it| если:
             *  1) элемент |edge_it| добавлен в |crossed_edges| не в текущей
             *  точке-событии и добавлено нижнее ребро
             *  2) группа не пустая и удалено нижнее ребро
             */

            if( (!new_edge_ovlp && (change_flags & LOWER_ADDED)) ||
                    (!empty_edge_ovlp && (change_flags & LOWER_DELETED)) )
            {
                if(edge_ovlp->inner_interval_it != nil)
                        throw "append_inner_intervals(): attempt to append "
                                "one edge twice";

                edge_ovlp->inner_interval_it =
                        inner_intervals.append(edge_it);
                last_interval = edge_it;
            }

            // Удаляем группу из |crossed_edges|, если она пуста
            if(empty_edge_ovlp)
            {
                if(edge_ovlp->inner_interval_it != nil)
                    throw "append_inner_intervals(): attempt to delete edge"
                            " present in |inner_intervals|";

                seq_item it = crossed_edges.succ(edge_it);
                crossed_edges.del_item(edge_it);
                edge_it = it;
            }
            else
            {
                // Сбрасываем флаги произведенных изменений
                edge_ovlp->flags &= ~CHANGES_MASK;
                edge_it = crossed_edges.succ(edge_it);
            }
        }
        else
            edge_it = crossed_edges.succ(edge_it);
    }
}

/******************************************************************************
 * Вычисляет возможные пересечения между |edge_it_pred| и следующим за ним ребром,
 * и между |edge_it_succ| и предшествующим ему ребром
 */
void y_trapez_sweep_box::compute_intersections()
{
    // Если |edge_it_pred| до сих пор |NIL|, значит нет ребер заканчивающихся,
    // проходящих или начинающихся в этой точке (т.е. |pnt_event|
    // изолированная точка); в этом случае больше ничего делать не надо.
    // Иначе нужно вычислить возможные пересечение между новыми соседями
    // в упорядоченном списке |crossed_edges|

    if (edge_it_pred != nil)
    {
        // У ребра |edge_it_pred| поменялось ребро, следующее за ним
        // в упорядоченном списке |crossed_edges|.
        // Обнуляем указатель на точку пересечения со следующим
        // за |edge_it_pred| ребром и вычисляем возможные пересечения.

        if (!use_optimization)
        {
            crossed_edges.change_inf(edge_it_pred, seq_item(nil));
            compute_intersection(event_points, crossed_edges, edge_it_pred);
            edge_it = crossed_edges.pred(edge_it_succ);
            if (edge_it != edge_it_pred)
                compute_intersection(event_points, crossed_edges, edge_it);
        }
        else
        {
            seq_item xit = crossed_edges.inf(edge_it_pred);
            if (xit)
            {
                edge_overlap_h edge_ovlp1 = crossed_edges.key(edge_it_pred);
                edge_overlap_h edge_ovlp2 = crossed_edges.key(edge_it_pred_succ); // edge_it_first
                inter_dic(edge_ovlp1, edge_ovlp2) = xit;
                crossed_edges.change_inf(edge_it_pred, seq_item(nil));
            }

            compute_intersection(event_points, crossed_edges, inter_dic, edge_it_pred);
            edge_it = crossed_edges.pred(edge_it_succ);
            if (edge_it != edge_it_pred)
                compute_intersection(event_points, crossed_edges, inter_dic, edge_it);
        }
    }
}

// Выполняет согласованное обновление счетчиков |inner_cnt|
// в связи с удалениями/добавлениями ребер в точках с одинаковым X-ми
void y_trapez_sweep_box::update_inner_counters(double sweep_x)
{
    // Цикл по всем интервалам обновления счетчиков |inner_cnt|
    while(!inner_intervals.empty())
    {
        seq_item start_edge_ovlp_it = inner_intervals.front();      // Начало очередного интевала

        // Получаем две соседние группы наложившихся ребер
        seq_item edge_ovlp_it0 = start_edge_ovlp_it;
        seq_item edge_ovlp_it = crossed_edges.succ(edge_ovlp_it0);
        edge_overlap_h edge_ovlp0 = crossed_edges.key(edge_ovlp_it0);
        edge_overlap_h edge_ovlp = crossed_edges.key(edge_ovlp_it);

        // Выполняем цикл пока значения счетчиков
        // |edge_ovlp->inner_cnt| и |edge_ovlp0->inner_cnt| не согласуются
        while(true)
        {
            if(edge_ovlp_it0 == start_edge_ovlp_it)
            {
                // Удаляем элемент из списка |inner_intervals|
                edge_ovlp0->inner_interval_it = nil;
                inner_intervals.del_item(inner_intervals.first_item());

                // Получаем следующий интервал
                start_edge_ovlp_it = (!inner_intervals.empty() ? inner_intervals.front() : nil);
            }

            // Вычисляем какое должно быть значение счетчика |inner_cnt| у |edge_ovlp|
            assert(0 <= edge_ovlp0->inner_cnt && 0 <= edge_ovlp0->edge_cnts[LOWER]);
            assert(0 <= edge_ovlp->edge_cnts[UPPER]);
            int due_inner_cnt = edge_ovlp0->inner_cnt +
                    edge_ovlp0->edge_cnts[LOWER] - edge_ovlp->edge_cnts[UPPER];

            if(edge_ovlp->inner_cnt != due_inner_cnt)
            {
                edge_ovlp->inner_cnt = due_inner_cnt;

                // Устанавливаем значения флагов "граница объединения"
                set_uborder_flags(edge_ovlp, sweep_x);

                // Продвигаемся дальше
                edge_ovlp_it0 = edge_ovlp_it;
                edge_ovlp_it = crossed_edges.succ(edge_ovlp_it);
                edge_ovlp0 = edge_ovlp;
                edge_ovlp = crossed_edges.key(edge_ovlp_it);
            }
            else
                break;
        }
    }
}

/***************************************************************
 * Вычисляет невертикальные отрезки границы объединения трапеций
 */
void y_trapez_sweep_box::unite_y_trapezoids(
    const list<y_trapezoid> &trapezoids,
    list<segment> &uborder_segms
    )
{
    // Инициализируем структуры данных
    init(trapezoids);
    p_union_edges = &uborder_segms;

    // Текущее положение заметающей линии
    double sweep_x = pnt_event.xcoord();

    // Основой цикл
    while (!event_points.empty())
    {
        // Переходим к следущей точке-событию
        event_it = event_points.min();
        pnt_event = event_points.key(event_it);
        cmp.set_event(pnt_event);

        if( !equal_eps(pnt_event.xcoord(),  sweep_x) )
        {
            // Новое положение заметающей линии
            // Выполняем согласованное обновление счетчиков |inner_cnt|
            // в связи с удалениями/добавлениями ребер в точках
            // с предыдущим X-ом
            update_inner_counters(sweep_x);

            sweep_x = pnt_event.xcoord();
        }

        // Получаем один из сегментов, проходящих через эту точку
        edge_it = event_points.inf(event_it);
        if (edge_it == nil)
        {
            // Нет ни точек пересечений, ни заканчивающихся ребер,
            // может быть только одна группа наложившихся проходящих ребро,
            // на которую попало начало нового ребра.
            // Если так, то находим её, чтобы корректно определить
            // |edge_it_succ| и |edge_it_pred|
            tmp_edge_ovlp->s = segment(pnt_event, pnt_event);
            edge_it = crossed_edges.lookup(tmp_edge_ovlp);
        }

        // Обнуляем указатели на ребро, следующее за пересекающимися
        // в текущей точке, и предшествующее им
        edge_it_succ = edge_it_pred = edge_it_pred_succ = nil;

        // Значение |NIL| для |edge_it_succ| и |edge_it_pred| после этих вычислений
        // означает, что нет сегментов оканчивающихся или проходящих через |pnt_event|

        if (edge_it != nil)     // key(edge_it) это заканчивающийся или проходящий сегмент
        {
            /*
             * Сначала нужно из счетчиков вычесть наложения удаляемых ребер
             * правых треугольников.
             * Затем нужно пересчитать счетчики наложений внутренностями
             * в связи с пересечением ребер.
             * После этого нужно из счетчиков вычесть наложения остальных
             * удаляемых ребер.
             */

            // Вычитаем из счетчиков наложения ребра правых треугольников,
            // заканчивающиеся в текущей точке. Остальные удаляемые ребра,
            // заносим в список |event_del_edges|
            prepare_del_edges();

            // Пересчитываем счетчики |inner_cnt| у ребер, пересекающихся в текущей
            // точке-событии, меняем порядок следования пересекшихся ребер на обратный.
            // Определяем |edge_it_succ| и |edge_it_pred|
            intersect_edges();

            // Вычитатаем из счетчиков наложений остальные удаляемые ребра
            delete_edges();
        }

        // Вставляем ребра, начинающиеся в |pnt_event|,
        // в список пересекаемых заметающей прямой
        insert_new_edges();

        // На флагов произведенных изменений в |edge_overlap_rep::flags|
        // формируем интервалы ребер в последовательности |crossed_edges|,
        // на которых либо все счетчики |inner_cnt| изменились,
        // либо все не изменились.
        // Удаляем из списка |crossed_edges| пустые группы наложившихся ребер
        append_inner_intervals();

        // В каждой точке-событии меняется ребро следующее за |edge_it_pred|
        // и ребро предшествующее |edge_it_succ|
        // Вычисляем возможные пересечения между новыми соседями
        compute_intersections();

        // Удаляем событие
        event_points.del_item(event_it);
    }

    assert(crossed_edges.size() == 2);  // Остались только ограничитель
    assert(edge_queue.size() == 1 && edge_del_queue.size() == 1);
}

/* Приводит трапецию с основаниями, параллельными оси Y в нормальный вид */
void y_trapezoid::normalize()
{
    // Убеждаемся что ребра невертикальные
    if(  ( equal_eps(edges[LOWER].dx(), 0.0) && !equal_eps(edges[LOWER].dy(), 0.0) )
            || ( equal_eps(edges[UPPER].dx(), 0.0) && !equal_eps(edges[UPPER].dy(), 0.0) )  )
        throw "y_trapezoid::normalize(): vertical edge";

    // Убеждаемся что ребра направлены слева направо
    if( less_eps(edges[LOWER].dx(), 0.0) )
       edges[LOWER] = segment(edges[LOWER].target(), edges[LOWER].source());

    if( less_eps(edges[UPPER].dx(), 0.0) )
       edges[UPPER] = segment(edges[UPPER].target(), edges[UPPER].source());

    // Проверяем что основания параллельны оси Y
    if(  !equal_eps( edges[LOWER].xcoord1(), edges[UPPER].xcoord1() ) ||
            !equal_eps( edges[LOWER].xcoord2(), edges[UPPER].xcoord2() )  )
        throw "y_trapezoid::normalize(): base of trapezoid is not parallel to axis Y";

    // Проверяем взаимное расположение нижней и верхней сторон
    if(  !less_eq_eps( edges[LOWER].ycoord1(), edges[UPPER].ycoord1() )||
            !less_eq_eps( edges[LOWER].ycoord2(), edges[UPPER].ycoord2() )  )
        throw "y_trapezoid::normalize(): lower and upper edges of trapezoid are permuted";
}

edge_overlap_rep::edge_overlap_rep()
{
    inner_cnt = INT_MIN;    // Чтобы отличать от отриц. из-за ошибки значений
    edge_cnts[LOWER] = edge_cnts[UPPER] = 0;
    flags = 0;
    inner_interval_it = nil;
}

edge_overlap_rep::edge_overlap_rep(const segment &s)
{
    this->s = s;
    inner_cnt = INT_MIN;    // Чтобы отличать от отриц. из-за ошибки значений
    edge_cnts[LOWER] = edge_cnts[UPPER] = 0;
    flags = 0;
    inner_interval_it = nil;
}

int compare(const edge_overlap_h &, const edge_overlap_h &)
{
    throw "compare(const edge_overlap_h &, const edge_overlap_h &): use 'edge_overlap_cmp(...)' instead";
    return 0;
}

/* Вычисляет невертикальные отрезки границы объединения трапеций */
void unite_y_trapezoids(const list<y_trapezoid> &trapezoids, list<segment> &uborder_segms)
{
    y_trapez_sweep_box sweep_box;
    sweep_box.unite_y_trapezoids(trapezoids, uborder_segms);
}

// Вынес определение, чтобы функция не была inline, чтобы её можно было вызывать при отладке
edge_overlap_rep* edge_overlap_h::ptr() const { return (edge_overlap_rep*)PTR; }