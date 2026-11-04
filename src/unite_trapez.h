#ifndef _TRAPEZ_UNION_H
#define	_TRAPEZ_UNION_H

// Флаги участка границы объединения
#define LOWER_BORDER   LOWER_EDGE
#define UPPER_BORDER   UPPER_EDGE
#define UNION_BORDER    2
#define BORDER_MASK    (1+2)

// Флаги произведенных изменений в группе наложившихся ребер
#define LOWER_DELETED   4
#define LOWER_ADDED     8
#define UPPER_DELETED   16
#define UPPER_ADDED     32
#define CHANGES_MASK    (4+8+16+32)

/************************************************************
 * Набор накладывающихся друг на друга верхних и нижних ребер
 * трапеций с основаниями параллельными оси Y
 */
struct edge_overlap_rep : public handle_rep
{
    // Наименьший сегмент, содержащий все накладывающиеся ребра
    segment s;

    // Кол-во трапеций, накрывающих участок сегмента своей внутренностью
    int inner_cnt;

    // Кол-во трапеций, накладывающихся на участок своим нижним/верхним ребром
    int edge_cnts[2];

    /* Флаги:
     *  UPPER_UBORDER – является верхней границей объединения трапеций
     *  LOWER_UBORDER – является нижней границей объединения трапеций
     *  UBORDER – является (верхней или нижней) границей объединения трапеций
     *  (на каждом интервале Y-ков участок сегмента может быть либо верхней границей,
     *  либо нижней, но не обеими одновременно)
     */
    int flags;

    // Начиная с какой точки является (следующим) участком границы
    point uborder_since;

    list_item inner_interval_it;    // Позиция группы наложившихся ребер в списке
                                    // |y_trapez_sweep_box::inner_intervals|
                                    // (нужна для диагностики проблем)

    edge_overlap_rep();

    edge_overlap_rep(const segment &s);
    
    ~edge_overlap_rep() {};
};

/*******************************************************************************
 * Всегда не равный NULL указатель на |edge_overlap_rep|, ведущий подсчет ссылок
 */
class edge_overlap_h : public handle_base
{
public:

    IMPLEMENT_HANDLE(edge_overlap_h, edge_overlap_rep)

    edge_overlap_h(const segment &s) {PTR = new edge_overlap_rep(s);}
};


/****************************************************
 *  Элемент упорядоченной последовательнсти с флагами
 */
struct flag_seq_item
{
    seq_item it;
    BYTE flags;

    flag_seq_item()
    {it = nil; flags = 0;}

    flag_seq_item(seq_item it, BYTE flags)
    {this->it = it; this->flags = flags;}
};

#endif	/* _TRAPEZ_UNION_H */
