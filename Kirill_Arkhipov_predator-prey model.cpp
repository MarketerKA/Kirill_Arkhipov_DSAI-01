//Kirill Arkhipov
//DSAI-01
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

class Victim {
public:
    int number_of_victims;
    double a1;
    double b1;
};

class Killers {
public:
    int number_of_killers;
    double a2;
    double b2;
};

double
victums_with_time(double victims, double time, double killers, double a1,
                  double b1, double a2, double b2) {
    return (victims - a2 / b2) * cos(sqrt(a2 * a1) * time) - (killers - a1 / b1) *
                                                                    (sqrt(a2) * b1 / (b2 * sqrt(a1))) *
                                                                    sin(sqrt(a1 * a2) * time) + a2 / b2;
}

double kill_with_time(double victims, double time, double killers, double a1,
                      double b1, double a2, double b2) {
    return (victims - a2 / b2) * (sqrt(a1) * b2) / (b1 * sqrt(a2)) * sin(sqrt(a1 * a2) * time) +
           (killers - a1 / b1) * cos(sqrt(a1 * a2) * time) + a1 / b1;
}

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

int main() {
#ifdef WIN32
    FILE *pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif
    Victim victim;
    Killers killers;
    double time;
    int num_points;
    cin >> victim.number_of_victims >> killers.number_of_killers;
    cin >> victim.a1 >> victim.b1;
    cin >> killers.a2 >> killers.b2;
    cin >> time >> num_points;

    fprintf(pipe, "%s\n", "plot '-' title 'Data' with lines, '-' title 'Predator-prey model' with lines");

    double coefficient = time / num_points;
    double time_work = 0.00;
    vector<double> v, k, t;
    while (time_work <= time) {
        double kil = kill_with_time(victim.number_of_victims, time_work,
                                    killers.number_of_killers, victim.a1, victim.b1, killers.a2, killers.b2);
        k.push_back(kil);
        t.push_back(time_work);
        double vic = victums_with_time(victim.number_of_victims, time_work,
                                       killers.number_of_killers, victim.a1, victim.b1, killers.a2, killers.b2);
        v.push_back(vic);


        time_work += coefficient;
    }
    cout << "t:\n";
    for (double l: t) {
        cout << fixed << setprecision(2) << l << " ";
    }
    cout << "\n" << "v:\n";
    int count = 0;
    for (double l: v) {
        fprintf(pipe, "%f\t%f\n", k[count], l);
        count++;
        cout << l << fixed << setprecision(2) << " ";
    }
    cout << "\n" << "k:\n";

    fprintf(pipe, "%s\n", "e");
    for (double l: k) {
        cout << fixed << setprecision(2) << l << " ";
    }
    fflush(pipe);
    pclose(pipe);

    return 0;
}