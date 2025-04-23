//
// Created by Riva Alkahal on 16/04/2025.
//
#include <iostream>
#include <unordered_map>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <Eigen/Dense>

using MatrixMap = std::unordered_map<std::string, std::vector<Eigen::MatrixXd>>;

std::mutex mtx;
std::condition_variable cv;
MatrixMap matrix_data;
bool done = false;

void computeMergedColumns(const MatrixMap& snapshot, int last_processed) {
    auto it = snapshot.find("design");
    if (it == snapshot.end()) return;

    const auto& matrices = it ->second;
    int upTo = matrices.size();  // process all available

    std::vector<int> cols = {0, 2};  // columns to extract
    int rows = matrices[0].rows();
    int totalCols = upTo * cols.size();

    Eigen::MatrixXd merged(rows, totalCols);
    int colPos = 0;

    for (int i = 0; i < upTo; ++i) {
        for (int j : cols) {
            merged.col(colPos++) = matrices[i].col(j);
        }
    }

    std::cout << "[Consumer] Merged size after " << upTo << " matrices: "
              << merged.rows() << "x" << merged.cols() << std::endl;
}

void producer(int total) {
    for (int i = 0; i < total; ++i) {
        Eigen::MatrixXd mat = Eigen::MatrixXd::Constant(3, 4, i);

        {
            std::lock_guard<std::mutex> lock(mtx);
            matrix_data["design"].push_back(mat);
        }

        cv.notify_one();
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }

    {
        std::lock_guard<std::mutex> lock(mtx);
        done = true;
    }

    cv.notify_all();
}

void consumer() {
    int last_processed = 0;

    while (true) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]() {
            return done || (matrix_data["design"].size() >= last_processed + 100);
        });

        if (matrix_data["design"].size() >= last_processed + 100) {
            MatrixMap snapshot = matrix_data;
            last_processed += 100;
            lock.unlock();

            computeMergedColumns(snapshot, last_processed);
        } else if (done) {
            break;
        }
    }

    std::cout << "[Consumer] Done.\n";
}

int main() {
    std::thread prod(producer, 500);
    std::thread cons(consumer);

    prod.join();
    cons.join();

    std::cout << "Main finished.\n";
    return 0;
}
