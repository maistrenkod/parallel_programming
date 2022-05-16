#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#define STRUCT_T_PEC 234
#define EXPECT_NUM_ARG 3
#define EXPECT_ST_INT 1000
#define START 0.001
#define END 1.0
#define NUM_TREND_INT 319
#define PI 3.14159265358979323846 



static inline double func(double x)
{
    return sin(1 / x);
};

static inline void print_error(const char* message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
};

long get_integer_number_from_arguments(char** argv, unsigned int pos)
{
    const int base = 10;
    char* endptr, * str;
    long val;

    str = argv[pos];

    errno = 0;
    val = strtol(str, &endptr, base);

    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0))
    {
        print_error("Bad strtol");
    }

    if (endptr == str)
        print_error("No digits were found in N");

    if (*endptr != '\0')
        print_error("Further characters after N");

    return val;
};

double get_double_number_from_arguments(char** argv, unsigned int pos)
{
    char* endptr, * str;
    double val;

    str = argv[pos];

    errno = 0;
    val = strtod(str, &endptr);

    if (endptr == str)
        print_error("No digits were found in N");

    if (*endptr != '\0')
        print_error("Further characters after N");

    return val;
};

double get_precision_from_args(int argc, char** argv)
{
    if (argc != EXPECT_NUM_ARG)
        print_error("Usage [num_procs] [precision]");

    const double val = get_double_number_from_arguments(argv, 2);

    if (val <= 0)
        print_error("precision should be a positive number");

    return val;
};

int get_num_procs_from_args(int argc, char** argv)
{
    if (argc != EXPECT_NUM_ARG)
        print_error("Usage [num_procs] [precision]");

    const long val = get_integer_number_from_arguments(argv, 1);

    if (val <= 0)
        print_error("num_procs should be a positive number");

    return (int)val;
};


double global_sum = 0;
pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;


typedef struct thread_data_t
{
    double stan_gran_step;
    size_t stan_gran_start;
    size_t stan_gran_num_steps;
    size_t dyn_gran_start;
    size_t dyn_gran_num_steps;
    double dyn_gran_start_point;
    double dyn_gran_end_point;
    int last_stan_step;
} thread_data_t;


void* calculate(void* args)
{
    const thread_data_t data = *(thread_data_t*)args;

    const double stan_gran_step = data.stan_gran_step;
    const size_t stan_gran_start = data.stan_gran_start;
    const size_t stan_gran_num_steps = data.stan_gran_num_steps;

    const size_t dyn_gran_start = data.dyn_gran_start;
    const size_t dyn_gran_num_steps = data.dyn_gran_num_steps;

    const double dyn_gran_end_point = data.dyn_gran_end_point;

    const int last_stan_step = data.last_stan_step;

    double sum = 0;
    double current_point = 0;
    double current_step = 0;

    if (dyn_gran_start == 0){
        current_step = 1 / (PI / 2 + (NUM_TREND_INT - 2) * PI) - START;
    }
    else{
        current_step = 1 / (PI / 2 + (NUM_TREND_INT - 2 - dyn_gran_start) * PI) - 1 / (PI / 2 + (NUM_TREND_INT - 1 - dyn_gran_start) * PI);
    }


    if (dyn_gran_start == 0){
        current_point = START;
    }
    else{
        current_point = 1 / (PI / 2 + (NUM_TREND_INT - 1 - dyn_gran_start) * PI);
    }

    size_t i = 0;
    for (i = dyn_gran_start; i < dyn_gran_num_steps + dyn_gran_start; i++){
        sum += (func(current_point) + 0.5 * current_step * func(current_point + current_step));
        current_point = 1 / (PI / 2 + (NUM_TREND_INT - 1 - i) * PI);
        current_step = 1 / (PI / 2 + (NUM_TREND_INT - 2 - i) * PI) - 1 / (PI / 2 + PI * (NUM_TREND_INT - 1 - i));
    }

    const size_t standard_granulation_end = (last_stan_step) ? (stan_gran_start + stan_gran_num_steps - 1) : (stan_gran_start + stan_gran_num_steps);
    current_point = dyn_gran_end_point + (stan_gran_start)*stan_gran_step;

    for (i = stan_gran_start; i < standard_granulation_end; i++){
        sum += (func(current_point) + func(current_point + stan_gran_step)) * 0.5 * stan_gran_step;
        current_point += stan_gran_step;
    }


    if (last_stan_step){
        const double last_standart_step = END - current_point;
        sum += (func(current_point) + func(current_point + last_standart_step)) * 0.5 * last_standart_step;
    }

    pthread_mutex_lock(&global_mutex);
    global_sum += sum;
    pthread_mutex_unlock(&global_mutex);

    return NULL;
}



int main(int argc, char** argv)
{
    const int num_procs = get_num_procs_from_args(argc, argv);
    const double precision = get_precision_from_args(argc, argv);

    pthread_t* pthreads = (pthread_t*)calloc(sizeof(*pthreads), num_procs);
    if (pthreads == NULL)
        print_error("Bad alloc");

    thread_data_t* pthread_data = (thread_data_t*)calloc(sizeof(*pthread_data), num_procs);
    if (pthread_data == NULL)
        print_error("Bad alloc");

    const double expected_step = sqrt(12.0 / EXPECT_ST_INT * precision / (END - START));

    double current_step = 1 / (PI / 2 + (NUM_TREND_INT - 2) * PI) - START;

    size_t dyn_gran_num_steps = 0;

    double standard_dynamic_step_threshold = 2 * expected_step;

    while (current_step < standard_dynamic_step_threshold && dyn_gran_num_steps < NUM_TREND_INT - 1){
        dyn_gran_num_steps++;
        current_step = 1 / (PI / 2 + (NUM_TREND_INT - 2 - dyn_gran_num_steps) * PI) - 1 / (PI / 2 + (NUM_TREND_INT - 1 - dyn_gran_num_steps) * PI);
    }

    const double dyn_gran_end_point = 1 / (PI / 2 + (NUM_TREND_INT - 2 - dyn_gran_num_steps) * PI);

    size_t stan_gran_num_steps = 0;
    if (dyn_gran_num_steps != 0){
        stan_gran_num_steps = (size_t)((END - 1 / (PI / 2 + (NUM_TREND_INT - 1 - dyn_gran_num_steps) * PI)) / expected_step);
    }
    else{
        stan_gran_num_steps = (size_t)((END - START) / expected_step);
    }

    size_t global_standard_granulation_distrib_counter = stan_gran_num_steps + 1;
    size_t global_dynamic_granulation_distrib_counter = dyn_gran_num_steps;

    const size_t ordinary_granulation_per_thread = (global_standard_granulation_distrib_counter + global_dynamic_granulation_distrib_counter) / (size_t)num_procs;

    int i = 0;
    for (i = 0; i < num_procs; i++){
        size_t standard_distrib_counter = 0;
        size_t dynamic_distrib_counter = 0;

        size_t dynamic_start = 0;
        size_t standard_start = 0;

        int last_stan_step = 0;

        const size_t target_distrib_counter_value = (i != num_procs - 1) ? ordinary_granulation_per_thread : (stan_gran_num_steps + 1 + dyn_gran_num_steps - ordinary_granulation_per_thread * (num_procs - 1));

        //dynamic granulation
        if (global_dynamic_granulation_distrib_counter != 0)
        {
            dynamic_start = (dyn_gran_num_steps - global_dynamic_granulation_distrib_counter);

            if (target_distrib_counter_value > global_dynamic_granulation_distrib_counter)
            {
                dynamic_distrib_counter += global_dynamic_granulation_distrib_counter;
                global_dynamic_granulation_distrib_counter = 0;
            }
            else
            {
                dynamic_distrib_counter += target_distrib_counter_value;
                global_dynamic_granulation_distrib_counter -= target_distrib_counter_value;
            }
        }

        // standard granulation
        if (dynamic_distrib_counter != target_distrib_counter_value){
            standard_start = (stan_gran_num_steps + 1 - global_standard_granulation_distrib_counter);

            standard_distrib_counter += (target_distrib_counter_value - dynamic_distrib_counter);
            global_standard_granulation_distrib_counter -= standard_distrib_counter;
        }

        if (global_standard_granulation_distrib_counter == 0){
            last_stan_step = 1;
        }

        pthread_data[i] = (thread_data_t){
            .stan_gran_step = expected_step,
            .stan_gran_start = standard_start,
            .stan_gran_num_steps = standard_distrib_counter,

            .dyn_gran_start = dynamic_start,
            .dyn_gran_num_steps = dynamic_distrib_counter,

            .dyn_gran_end_point = dyn_gran_end_point,

            .last_stan_step = last_stan_step
        };

        if (pthread_create(&pthreads[i], NULL, calculate, &pthread_data[i]) == -1) {
            print_error("Bad pthread create");
        }

    }

    for (i = 0; i < num_procs; i++){
        if (pthread_join(pthreads[i], NULL) == -1)
            print_error("Bad pthread join");
    }

    printf("Integral: %lg\n", global_sum);

    return 0;
}