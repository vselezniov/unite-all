#include <cstdlib>
#include <ctime>
#include <cassert>
#include <iostream>

using namespace std;

#include "LEDA/point.h"
#include "LEDA/segment.h"
#include "LEDA/gen_polygon.h"
#include "LEDA/dictionary.h"

#include "float_compare.h"
#include "unite_all.h"

const char USAGE[] =    "Unites regions formed by N repetitions of simple polygon\n"
                        "Using:\n"
                        "\t unite-all N\n";

#define dimof(x)    (sizeof(x)/sizeof(x[0]))

// Если |UNITE_ALL| не определено, то используется LEDA реализация
#define UNITE_ALL

// Удаляет последовательные вершины, лежащие на одной прямой
gen_polygon eliminate_colinear_vertices(const gen_polygon &gpoly);

// Вычисляет объединение вереницы простых частично накладывающихся многоугольников
// (растянутой вдоль оси Y)
int unite_overlapping_train_y(int argc, char** argv);

// Вычисляет объединение сцепленных "расчесок" (растянутых вдоль оси Y)
int unite_combs_y(int argc, char** argv);

// Вычисляет объединение многих простых многоугольников (растянутых вдоль оси Y)
int unite_many_simple_y(int argc, char** argv);

// Вычисляет объединение многоугольных регионов (растянутых вдоль оси Y)
int unite_regions_y(int argc, char** argv);

// Выполняет простой тест объединения многоугольников
int unite_simple_test();

int main(int argc, char** argv)
{
    return unite_overlapping_train_y(argc, argv);

    //return unite_combs_y(argc, argv);

    //return unite_many_simple_y(argc, argv);

    //return unite_regions_y(argc, argv);

    //return unite_simple_test();

    return 0;
}

// Вычисляет объединение вереницы простых частично накладывающихся многоугольников
// (растянутой вдоль оси Y)
int unite_overlapping_train_y(int argc, char** argv)
{
    int N;          // Кол-во повторений простых многоугольников

    point point_array[] = { point(-1, -3), point( 1, -3), point( 3, -1), point( 3,  1),
                            point( 1,  3), point(-1,  3), point(-3,  1), point(-3, -1)};

    list<point> point_list;                 // Точки простого многоугольника
    gen_polygon poly;
    list<gen_polygon> poly_list;

    gen_polygon res_poly;

    clock_t clockStart;
    double dUniteSecs;
    int nUnionVertN;

    if( argc != 2 || !(N = atoi(argv[1])) )
    {
        cout << USAGE;
        goto err_lbl;
    }

    cout << "overlapping_train" << endl;

    point_list = list<point>(point_array, point_array + dimof(point_array));

    poly = gen_polygon(point_list).translate(3, 3);
    for(int i = 0; i < N; ++i)
    {
        if(0 < i)
            poly = poly.translate(0, 4);
        poly_list.append(poly);
    }

    // Объединяем многоугольники
    clockStart = clock();

#ifdef  UNITE_ALL
    res_poly = unite_all(poly_list);
#else
    forall(poly, poly_list)
        res_poly = res_poly.unite(poly);
#endif

    dUniteSecs = double(clock() - clockStart) / CLOCKS_PER_SEC;

    //cout << "\nInput polygons:\n" << poly_list;
    //cout << "\nReslut of uniting:\n" << res_poly << endl;

    assert(res_poly.polygons().size() == 1);
    nUnionVertN = res_poly.polygons().front().vertices().size();

    cout << "union vertex number: " << nUnionVertN << (nUnionVertN == 6*N + 2 ? "" : "wrong") << endl;
    cout << "time: uniting " << dUniteSecs << " sec." << endl;

    return 0;

err_lbl:
    return -1;
}

// Вычисляет объединение сцепленных "расчесок" (растянутых вдоль оси Y)
int unite_combs_y(int argc, char** argv)
{
    int N;                      // Кол-во зубьев у расчески

    list<point> comb1_points, comb2_points;
    gen_polygon comb1, comb2;
    list<gen_polygon> poly_list;

    gen_polygon res_poly;
    const list<polygon> *p_res_contours;

    clock_t clockStart;
    double dUniteSecs;
    int nUnionVertN;

    list_item it;
    int i;

    if( argc != 2 || !(N = atoi(argv[1])) )
    {
        cout << USAGE;
        goto err_lbl;
    }

    cout << "combs" << endl;

    // Формируем первую "расческу"
    for(i = 0; i < N; ++i)
    {
        comb1_points.push_back(point(1, 1 + i*2));
        comb1_points.push_back(point(1, 2 + i*2));
        comb1_points.push_back(point(3, 2 + i*2));
        comb1_points.push_back(point(3, 3 + i*2));
    }
    --i;
    comb1_points.push_back(point(0, 3 + i*2));
    comb1_points.push_back(point(0, 0));
    comb1_points.push_back(point(3, 0));
    comb1_points.push_back(point(3, 1));

    comb1 = gen_polygon(comb1_points);
    //cout << comb1;

    // Формируем вторую "расческу"
    for(i = 0; i < N; ++i)
    {
        comb2_points.push_back(point(2, 1 + i*2));
        comb2_points.push_back(point(2, 2 + i*2));
        comb2_points.push_back(point(3, 2 + i*2));
        comb2_points.push_back(point(3, 3 + i*2));
    }
    --i;
    comb2_points.push_back(point(4, 3 + i*2));
    comb2_points.push_back(point(4, 0));
    comb2_points.push_back(point(3, 0));
    comb2_points.push_back(point(3, 1));
    comb2_points.reverse();

    comb2 = gen_polygon(comb2_points);
    //cout << comb2;

    poly_list.append(comb1);
    poly_list.append(comb2);

    // Объединяем "расчески"
    clockStart = clock();

#ifdef  UNITE_ALL
    res_poly = unite_all(poly_list);
#else
    res_poly = comb1.unite(comb2);
    res_poly = eliminate_colinear_vertices(res_poly);
#endif

    dUniteSecs = double(clock() - clockStart) / CLOCKS_PER_SEC;

    //cout << "\nInput polygons:\n" << poly_list;
    //cout << "\nReslut of uniting:\n" << res_poly << endl;

    p_res_contours = &res_poly.polygons();
    assert(p_res_contours->size() == N + 1);

    nUnionVertN = 0;
    for(i = 0, it = p_res_contours->first_item(); i <= N; ++i, it = p_res_contours->succ(it) )
    {
        const polygon &cont = p_res_contours->contents(it);
        assert(cont.vertices().size() == 4);
        assert(cont.orientation() == (i == 0 ? 1 : -1) );
        nUnionVertN += cont.size();
    }

    cout << "union vertex number: " << nUnionVertN << (nUnionVertN == 4*(N+1) ? "" : "wrong") << endl;
    cout << "time: " << dUniteSecs << " sec." << endl;

    return 0;

err_lbl:
    return -1;
}

// Вычисляет объединение многих простых многоугольников (растянутых вдоль оси Y)
int unite_many_simple_y(int argc, char** argv)
{
    int N;              // Кол-во повторений простых многоугольников

    point point_array[] = { point(-1, -3), point( 1, -3), point( 3, -1), point( 3,  1),
                            point( 1,  3), point(-1,  3), point(-3,  1), point(-3, -1)};

    list<point> point_list;              // Точки простого многоугольника
    gen_polygon poly;
    list<gen_polygon> poly_list;

    gen_polygon res_poly;

    clock_t clockStart;
    double dUniteSecs;
    int nUnionVertN;

    if( argc != 2 || !(N = atoi(argv[1])) )
    {
        cout << USAGE;
        goto err_lbl;
    }

    cout << "many_simple" << endl;

    point_list = list<point>(point_array, point_array + dimof(point_array));


    // Многоугольники 1-го региона
    poly = gen_polygon(point_list).translate(3, 6);
    for(int i = 0; i < N; ++i)
    {
        if(0 < i)
            poly = poly.translate(0, 7);
        poly_list.append(poly);
    }

    // Многоугольники 2-го региона
    poly = gen_polygon(point_list).translate(7, 3);
    for(int i = 0; i < N; ++i)
    {
        if(0 < i)
            poly = poly.translate(0, 7);
        poly_list.append(poly);
    }

    // Объединяем многоугольники
    clockStart = clock();

#ifdef  UNITE_ALL
    res_poly = unite_all(poly_list);
#else
    forall(poly, poly_list)
        res_poly = res_poly.unite(poly);
#endif

    dUniteSecs = double(clock() - clockStart) / CLOCKS_PER_SEC;

    //cout << "\nInput polygons:\n" << poly_list;
    //cout << "\nReslut of uniting:\n" << res_poly << endl;

    assert(res_poly.polygons().size() == 1);
    nUnionVertN = res_poly.polygons().front().vertices().size();

    cout << "union vertex number: " << nUnionVertN << (nUnionVertN == 12*N + 2 ? "" : "wrong") << endl;
    cout << "time: " << "  uniting  " << dUniteSecs << " sec." << endl;

    return 0;

err_lbl:
    return -1;
}

// Вычисляет объединение многоугольных регионов (растянутых вдоль оси Y)
int unite_regions_y(int argc, char** argv)
{
    int N;          // Кол-во повторенний простых многоугольников

    point point_array[] = { point(-1, -3), point( 1, -3), point( 3, -1), point( 3,  1),
                            point( 1,  3), point(-1,  3), point(-3,  1), point(-3, -1)};

    list<point> point_list;                 // Точки простого многоугольника
    polygon poly;
    list<polygon> cont_list1, cont_list2;   // Контуры объединяемых регионов

    gen_polygon region1, region2;
    list<gen_polygon> poly_list;
    gen_polygon res_poly;

    clock_t clockStart;
    double dUniteSecs;
    int nUnionVertN;

    if( argc != 2 || !(N = atoi(argv[1])) )
    {
        cout << USAGE;
        goto err_lbl;
    }

    cout << "regions" << endl;

    point_list = list<point>(point_array, point_array + dimof(point_array));

    // Формируем контуры 1-го региона
    poly = polygon(point_list).translate(3, 6);
    for(int i = 0; i < N; ++i)
    {
        if(0 < i)
            poly = poly.translate(0, 7);
        cont_list1.append(poly);
    }

    // Формируем контуры 2-го региона
    poly = polygon(point_list).translate(7, 3);
    for(int i = 0; i < N; ++i)
    {
        if(0 < i)
            poly = poly.translate(0, 7);
        cont_list2.append(poly);
    }

    region1 = gen_polygon(cont_list1);
    region2 = gen_polygon(cont_list2);

    poly_list.append(region1);
    poly_list.append(region2);

    // Объединяем регионы
    clockStart = clock();
    
#ifdef  UNITE_ALL
    res_poly = unite_all(poly_list);
#else
    res_poly = region1.unite(region2);
#endif

    dUniteSecs = double(clock() - clockStart) / CLOCKS_PER_SEC;

    //cout << "\nInput polygons:\n" << poly_list;
    //cout << "\nReslut of uniting:\n" << res_poly << endl;

    assert(res_poly.polygons().size() == 1);
    nUnionVertN = res_poly.polygons().front().vertices().size();

    cout << "union vertex number: " << nUnionVertN << (nUnionVertN == 12*N + 2 ? "" : "wrong") << endl;
    cout << "time: " << "  uniting  " << dUniteSecs << " sec." << endl;

    return 0;

err_lbl:
    return -1;
}


// Выполняет простой тест объединения многоугольников
int unite_simple_test()
{
    {
        point point_array[5][4] = {
                                    { point(0, 0), point(6, 0), point(6, 1), point(0, 1) },
                                    { point(5, 0), point(6, 0), point(6, 5), point(5, 5) },
                                    { point(6, 4), point(6, 5), point(0, 5), point(0, 4) },
                                    { point(0, 5), point(0, 0), point(1, 0), point(1, 5) },
                                    { point(2, 2), point(4, 2), point(4, 3), point(2, 3) }
                                 };

        list<gen_polygon> poly_list;
        gen_polygon res_poly;

        for(int i = 0; i < 5; ++i)
        {
            list<point> point_list(point_array[i], point_array[i] + 4);
            poly_list.append(gen_polygon(point_list));
        }

#ifdef  UNITE_ALL
        res_poly = unite_all(poly_list);
#else   // LEDA
        gen_polygon poly;
        forall(poly, poly_list)
            res_poly = res_poly.unite(poly);
        res_poly = eliminate_colinear_vertices(res_poly);
#endif

        cout << "\nInput polygons:\n" << poly_list;
        cout << "\nReslut of uniting:\n" << res_poly << endl;
    }

    {
        point cont11[] = {  point(-7, -3), point(2, -3), point(2, 8), point(-7, 8) };
        point cont12[] = {  point(-5, 0), point(-5, 6), point(0, 6), point(0, 0) };

        point cont2[] = {   point(0, 0), point(-5, 0), point(-5, -6), point(7, -6), point(7, 4), point(-5, 4),
                            point(0, 0), point(0, 2), point(5, 2), point(5, -4), point(0, -4) };

        list<point> point_list;
        list<polygon> cont_list;
        list<gen_polygon> poly_list;
        gen_polygon res_poly;


        cont_list.clear();

        point_list.clear();
        point_list.insert(point_list.end(), cont11, cont11 + dimof(cont11));
        cont_list.append(polygon(point_list));

        point_list.clear();
        point_list.insert(point_list.end(), cont12, cont12 + dimof(cont12));
        cont_list.append(polygon(point_list));

        poly_list.append(gen_polygon(cont_list));


        cont_list.clear();

        point_list.clear();
        point_list.insert(point_list.end(), cont2, cont2 + dimof(cont2));
        cont_list.append(polygon(point_list));

        poly_list.append(gen_polygon(cont_list));


#ifdef  UNITE_ALL
        res_poly = unite_all(poly_list);
#else   // LEDA
        gen_polygon poly;
        forall(poly, poly_list)
            res_poly = res_poly.unite(poly);
        res_poly = eliminate_colinear_vertices(res_poly);
#endif

        cout << "\nInput polygons:\n" << poly_list;
        cout << "\nReslut of uniting:\n" << res_poly << endl;
    }

    return 0;
}

// Удаляет последовательные вершины, лежащие на одной прямой
gen_polygon eliminate_colinear_vertices(const gen_polygon &gpoly)
{
    list<polygon> contsRes;
    polygon contSrc;

    forall(contSrc, gpoly.polygons())
    {
        contsRes.push_back(contSrc.eliminate_colinear_vertices());
    }

    return gen_polygon(contsRes);
}