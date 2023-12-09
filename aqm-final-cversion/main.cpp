#include <iostream>
#include <utility>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <set>


class Lattice{
public:
    int N;
    std::vector<int> index = std::vector<int> {};
    std::vector<std::vector<int>> hop = std::vector<std::vector<int>> {};
    Lattice(){
        this->N = 0;
    }
    explicit Lattice(int N){
        this->N = N;
        this->index = std::vector<int> (N*N*N);
        this->hop.resize(7);
        for (auto i = 0; i < 7; ++i){
            this->hop[i].resize(N*N*N);
        }
        // Initialize the hopping matrix
        for (auto i = 0; i < std::pow(N, 3); ++i){
            this->index[i] = i;
            this->hop[0][i] = i;
            this->hop[1][i] = i + 1;
            this->hop[2][i] = i - 1;
            this->hop[3][i] = i + this->N;
            this->hop[4][i] = i - this->N;
            this->hop[5][i] = i + N*N;
            this->hop[6][i] = i - N*N;
        }
        // Fix the hopping matrix
        for (auto i = 0; i < N; ++i){
            for (auto j = 0; j < N; ++j){
                this->hop[1][j*N*N + i*N + N - 1] -= N;
                this->hop[2][j*N*N + i*N] += N;
                this->hop[3][j*N*N + (N-1) * N + i] -= (N*N);
                this->hop[4][j*N*N + i] += (N*N);
                this->hop[5][(N-1)*N*N + N*i + j] -= (N*N*N);
                this->hop[6][N*i + j] += (N*N*N);
            }
        }

    }
};


class SpaceTimePoint{
public:
    int time_index{};
    int space_index{};
    bool occupied{};
    SpaceTimePoint* head;
    SpaceTimePoint* tail;
    SpaceTimePoint()= default;
    SpaceTimePoint(int time_index_,
                   int space_index_,
                   SpaceTimePoint* head_,
                   SpaceTimePoint* tail_) {
        this->time_index = time_index_;
        this->space_index = space_index_;
        this->head = head_;
        this->tail = tail_;
        this->occupied = false;
    }
};

class GridWorld{
public:
    int len_time;
    int len_space;
    Lattice lattice;
    std::vector<std::vector<SpaceTimePoint*>> grid;
    GridWorld(){};
    GridWorld(int len_time_, Lattice lattice_){
        this->len_time = len_time_;
        this->lattice = std::move(lattice_);
        this->len_space = this->lattice.N*this->lattice.N*this->lattice.N;
        grid.resize(this->len_time);
        for (auto i = 0; i < this->len_time; ++i){
            grid[i].resize(this->len_space);
            for (auto j = 0; j < this->len_space; ++j){
                grid[i][j] = new SpaceTimePoint(i, j, nullptr, nullptr);
            }
        }
    }
    bool trajectory_initialization(){
        int space_index = rand() % this->len_space;
        int time_index = 0;
        auto node = this->grid[0][space_index];
        node->head = this->grid[this->len_time - 1][space_index];
        this->grid[this->len_time-1][space_index]->tail = node;
        node->occupied = true;

        for (auto i = 0; i < this->len_time; ++i){
            time_index = (time_index + 1) % (this->len_time);
            auto next_node = this->grid[time_index][space_index];
            next_node->occupied = true;
            next_node->head = node;
            node->tail = next_node;
            node = next_node;
        }
        return true;
    }
};

class Worm{
public:
    GridWorld gridWorld;
    int direction;
    int temporal_direction;
    SpaceTimePoint* position;
    double mu;
    double epsilon;
    double beta;
    bool keep_go = false;
    bool do_site_update = false;
    int forward_update = 0;

    Worm(){};
    Worm(GridWorld gridWorld_, int direction_, int temporal_direction_, SpaceTimePoint *position_, double mu_, double epsilon_, double beta_){
        gridWorld = std::move(gridWorld_);
        direction = direction_;
        temporal_direction = temporal_direction_;
        position = position_;
        mu = mu_;
        epsilon = epsilon_;
        beta = beta_;

    }
    void site_update(){
        if (this->temporal_direction == 1){
            this->forward_update += 1;
        }
        this->direction*= -1;
        double rand_num = (double) rand() / (double) RAND_MAX;
        if (rand_num > std::fmin(1, std::exp(this->temporal_direction * this->mu * this->epsilon))){
            this->temporal_direction *= -1;
        }
    }
    void backward_bond_update(){
        this->position->occupied = false;
        if(this->position->head == nullptr){
            this->keep_go = false;
            this->do_site_update = false;
            return;
        }
        auto new_position =this->position->head;
        this->position->head = nullptr;
        new_position->tail = nullptr;
        this->direction = 1;
        this->position = new_position;
        this->keep_go = true;
        this->do_site_update = true;
    }
    void forward_bond_update(){
        this->forward_update += 1;
        this->position->occupied = true;
        double rand_num = (double)rand() / (double) RAND_MAX;
        int hopping_direction = int(rand_num / this->epsilon + 1);
        if (hopping_direction > 6){
            hopping_direction = 0;
        }
        auto new_position = this->gridWorld.grid[
                (this->position->time_index + 1) % this->gridWorld.len_time
                ][this->gridWorld.lattice.hop[hopping_direction][this->position->space_index]
                  ];


        this->position->tail = new_position;
        if (new_position->head == nullptr){
            new_position->head = this->position;
            new_position->occupied = true;
            this->position = new_position;
            this->direction = 1;
            this->temporal_direction = 1;
            this->keep_go = true;
            this->do_site_update = true;
            return;
        }
        else {
            auto temp_position = new_position->head;
            temp_position->tail = nullptr;
            new_position->head = this->position;
            this->position = temp_position;
            new_position->occupied = true;
            this->direction = -1;
            this->temporal_direction = -1;
            this->keep_go = true;
            this->do_site_update = false;
            return;
        }
    }
    void update(){
        this->forward_update = 0;
        while(true){
            if (this->temporal_direction == 1) this->forward_bond_update();
            else this->backward_bond_update();
            if (this->do_site_update) this->site_update();
            if (not keep_go) return;
        }
    }
};

class WormAlgorithm{
public:
    Lattice lattice;
    double mu{};
    double epsilon{};
    double beta{};
    int time_length{};
    std::vector<double> e_list;
    std::vector<double> n_list;
    GridWorld grid_world;
    Worm worm;

    WormAlgorithm()= default;
    WormAlgorithm(Lattice lattice_, double mu_, double epsilon_, double beta_){
        this->lattice = std::move(lattice_);
        this->mu = mu_;
        this->epsilon = epsilon_;
        this->beta = beta_;
        this->time_length = int(beta / epsilon);
        this->e_list = std::vector<double> (time_length);
        this->n_list = std::vector<double> (time_length);
        this->grid_world = GridWorld(time_length, lattice);
    }
    bool worm_initialization(){
        int rand_t = rand() % this->time_length;
        int rand_r = rand() % this->grid_world.len_space;
        auto cursor = this->grid_world.grid[rand_t][rand_r];

        if(not cursor->occupied){
            double rand_num = (double) rand() / (double)RAND_MAX;
            if(rand_num < std::fmin(1, std::exp(this->mu * this->epsilon))){
                this->worm = Worm(this->grid_world, -1, 1, cursor, this->mu, this->epsilon, this->beta);
                return true;
            }
            else{
                return false;
            }
        }
        else{
            this->worm = Worm(this->grid_world, -1, -1, cursor, this->mu, this->epsilon, this->beta);
            return true;
        }
    }
    std::vector<double> calculate_expectation(){
        double n_spatial_hop = 0.;
        double n_temporal_hop = 0.;
        double n = 0.;
        for(auto space_index_start = 0 ; space_index_start < this->grid_world.len_space; ++ space_index_start){
            auto cursor = this->grid_world.grid[0][space_index_start];
            if(not cursor->occupied){
                continue;
            }
            n += 1;
            for (auto i = 0; i < this->time_length - 1; ++i){
                auto next = cursor->tail;
                if(next->space_index == cursor->space_index){
                    n_temporal_hop += 1;
                }
                else{
                    n_spatial_hop += 1;
                }
                cursor = next;
            }
        }
        auto k = this->beta / this->epsilon;
        auto e = (-1.0 / (k * this->epsilon)) * n_spatial_hop + (6. / (k * (1. - 6. * this->epsilon))) * n_temporal_hop;


        int winding_number = 0;
        auto checked = std::set<int> {};
        for(auto space_index_start = 0 ; space_index_start < this->grid_world.len_space; ++ space_index_start) {
            int one_winding = 0;
            auto cursor = this->grid_world.grid[0][space_index_start];
            if (not cursor->occupied) {
                continue;
            }
            if (checked.count(space_index_start)){
                continue;
            }

            auto start = cursor;
            int increment = 0;
            while(true){
                for (auto i = 0; i < this->time_length; ++i){
                    if(
                            (cursor->space_index / this->grid_world.lattice.N == 0) and
                            (cursor->tail->space_index % this->grid_world.lattice.N == (this->grid_world.lattice.N - 1))
                            ){
                            increment -= 1;
                        }
                    else if(
                            (cursor->space_index % this->grid_world.lattice.N == (this->grid_world.lattice.N - 1)) and
                            (cursor->tail->space_index % this->grid_world.lattice.N == 0)
                            ){
                        increment += 1;
                    }
                    else {
                        increment += (
                                (cursor->tail->space_index % this->grid_world.lattice.N)
                                - (cursor->space_index % this->grid_world.lattice.N)
                        );
                    }
                    cursor = cursor->tail;
                }
                if (cursor == start){
                    break;
                }
                else{
                    checked.insert(cursor->space_index);
                }
            }
//            std::cout << increment / this->lattice.N << ", " << increment % this->lattice.N << std::endl;
            winding_number += (increment / this->lattice.N);

        }
        return std::vector<double> {n, e, (double) winding_number, 0};
    }
    std::vector<double> single_run(){
        if(this->worm_initialization()){
            this->worm.update();
            auto expectation = this->calculate_expectation();
            expectation[3] = this->worm.forward_update;
            return expectation;
        }
        else{
            auto expectation = this->calculate_expectation();
            expectation[3] = 0.;
            return expectation;
        }
    }
    std::vector<std::vector<double>> run(int n_iter){
        this->grid_world.trajectory_initialization();
        auto result = std::vector<std::vector<double>>(n_iter);
        for (auto i = 0; i < n_iter; ++i){
            auto expectation_value = this->single_run();
            result[i]= expectation_value;
        }
        return result;
    }
};

int main(){
    int N = 8;
    auto lattice = Lattice(N);
    auto simulation = WormAlgorithm(lattice, 1., 0.01, 1.1);
    int n_iter = 30000;
    auto start = std::chrono::high_resolution_clock::now();

    auto result = simulation.run(n_iter);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "time: " << duration.count()  << "ms" << std::endl;
    std::ofstream result_file;
    result_file.open("result.txt");
    result_file << "N,E,W,X\n";
    for (auto & v : result){
        result_file << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "\n";
    }
    result_file.close();
}


