#ifndef __PROGTEST__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <array>
#include <iterator>
#include <set>
#include <list>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <compare>
#include <queue>
#include <stack>
#include <deque>
#include <memory>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <stdexcept>
#include <condition_variable>
#include <pthread.h>
#include <semaphore.h>
#include "progtest_solver.h"
#include "sample_tester.h"

using namespace std;

#include "atomic-logger.hpp"

#define LOG(LOGS) /*ATOMIC_LOG(tid << LOGS << std::endl)*/
#else
#define LOG(LOGS)
#endif /* __PROGTEST__ */

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef __PROGTEST__
#define SLEEP std::this_thread::sleep_for(std::chrono::milliseconds(rand()%10));
#else
#define SLEEP
#endif

template<typename T>
class CAtomicQueue {
public:
    void push(T value);

    T pop();

    queue<T> m_Queue;
    mutex m_Mutex;
    condition_variable m_ConditionNotEmpty;
};

template<typename T>
T CAtomicQueue<T>::pop() {
    unique_lock lock(m_Mutex);
    m_ConditionNotEmpty.wait(lock, [this] {
        return !m_Queue.empty();
    });
    auto front = std::move(m_Queue.front());
    m_Queue.pop();
    return front;
}

template<typename T>
void CAtomicQueue<T>::push(T value) {
    [[maybe_unused]] lock_guard lock(m_Mutex);
    m_Queue.push(std::move(value));
    m_ConditionNotEmpty.notify_all();
    SLEEP
}

class CEnhancedPack {
public:
    explicit CEnhancedPack(AProblemPack pack) : m_Pack(std::move(pack)) {}

    AProblemPack m_Pack;

    void waitSolved();

    void incrementSolved();

private:
    size_t m_SolvedCounter = 0;
    mutex m_Mutex;
    condition_variable m_CondVar;
};

class CSmartPolygon {
public:
    explicit CSmartPolygon(APolygon polygon, bool isCnt, shared_ptr<CEnhancedPack> enhancedPack)
        : m_Polygon(std::move(polygon)), m_IsCnt(isCnt), m_EnhancedPack(std::move(enhancedPack)) {}
    APolygon m_Polygon;
    bool m_IsCnt;
    shared_ptr<CEnhancedPack> m_EnhancedPack;
};

void CEnhancedPack::waitSolved() {
    unique_lock lock(m_Mutex);
    m_CondVar.wait(lock, [this]() {
        return m_SolvedCounter == m_Pack->m_ProblemsCnt.size() + m_Pack->m_ProblemsMin.size();
    });
    SLEEP
}

void CEnhancedPack::incrementSolved() {
    lock_guard lock(m_Mutex);
    m_SolvedCounter++;
    if (m_SolvedCounter == m_Pack->m_ProblemsCnt.size() + m_Pack->m_ProblemsMin.size()) m_CondVar.notify_all();
    SLEEP
}

mutex g_SolverMutex;

class CEnhancedSolver  {
public:
    CEnhancedSolver() = default;
    explicit CEnhancedSolver(AProgtestSolver solver, vector<shared_ptr<CEnhancedPack>> packs = {})
            : m_Solver(std::move(solver)),
              m_Packs(std::move(packs)) {}

    void solveAndMarkSolved();

    AProgtestSolver m_Solver;
    vector<shared_ptr<CEnhancedPack>> m_Packs;
};

void CEnhancedSolver::solveAndMarkSolved() {
    m_Solver->solve();
    for (auto &pack: m_Packs) {
        pack->incrementSolved();
    }
    SLEEP
}


class COptimizer {
public:
    static bool usingProgtestSolver() {
        return true;
    }

    static void checkAlgorithmMin(APolygon p) {
        // dummy implementation if usingProgtestSolver() returns true
    }

    static void checkAlgorithmCnt(APolygon p) {
        // dummy implementation if usingProgtestSolver() returns true
    }

    void start(int threadCount);

    void stop();

    void addCompany(ACompany company);

    void fetchNewPacks(const ACompany &company, const shared_ptr<CAtomicQueue<shared_ptr<CEnhancedPack>>> &ptr);

    void getSolvedPacks(const ACompany &company, const shared_ptr<CAtomicQueue<shared_ptr<CEnhancedPack>>> &ptr);

    void fillSolversAndSolve();

    CEnhancedSolver fetchSolver(AProgtestSolver (*createSolver)());

private:
    vector<thread> m_WorkingThreads;
    vector<thread> m_CompanyInputCommunicator;
    vector<thread> m_CompanyOutputCommunicator;
    CAtomicQueue<CSmartPolygon> m_PolygonsToSolve;
    atomic<int> m_CompanyCounter = 0;
    mutex m_Mutex;
    vector<ACompany> m_Companies;
    CEnhancedSolver m_MinSolver;
    CEnhancedSolver m_CntSolver;
};


void COptimizer::addCompany(ACompany company) {
    m_Companies.push_back(std::move(company));
    m_CompanyCounter++;
}

void COptimizer::start(int threadCount) {
    for (const auto &item: m_Companies) {
        shared_ptr queuePtr = make_shared<CAtomicQueue<shared_ptr<CEnhancedPack>>>();
        m_CompanyInputCommunicator.emplace_back(&COptimizer::fetchNewPacks, this, item, queuePtr);
        m_CompanyOutputCommunicator.emplace_back(&COptimizer::getSolvedPacks, this, item, std::move(queuePtr));
    }
    for (int i = 0; i < threadCount; ++i) {
        m_WorkingThreads.emplace_back(&COptimizer::fillSolversAndSolve, this);
    }
}

void COptimizer::stop() {
    for (auto &thr: m_CompanyInputCommunicator) {
        thr.join();
    }
    for (auto &thr: m_WorkingThreads) {
        thr.join();
    }
    for (auto &thr: m_CompanyOutputCommunicator) {
        thr.join();
    }
}

void COptimizer::fetchNewPacks(const ACompany &company, const shared_ptr<CAtomicQueue<shared_ptr<CEnhancedPack>>> &ptr) {
    while (true) {
        LOG("FETCH")
        auto pack = company->waitForPack();
        SLEEP
        if (pack == nullptr) {
            ptr->push(nullptr);
            m_PolygonsToSolve.push(CSmartPolygon(nullptr, false, nullptr));
            break;
        }
        auto package = make_shared<CEnhancedPack>(pack);
        ptr->push(package);
        lock_guard lock (m_PolygonsToSolve.m_Mutex);
        for (auto &polygon: pack->m_ProblemsCnt) {
            m_PolygonsToSolve.m_Queue.emplace(polygon, true, package);
        }
        for (auto &polygon: pack->m_ProblemsMin) {
            m_PolygonsToSolve.m_Queue.emplace(polygon, false, package);
        }
        SLEEP
        package = nullptr;
    }
}

void COptimizer::getSolvedPacks(const ACompany &company, const shared_ptr<CAtomicQueue<shared_ptr<CEnhancedPack>>> &ptr) {
    while (true) {
        LOG("GET")
        auto front = ptr->pop();
        if (front == nullptr) break;
        front->waitSolved();
        SLEEP
        company->solvedPack(front->m_Pack);
    }
}

void COptimizer::fillSolversAndSolve() {
    if (m_CompanyInputCommunicator.empty() && m_CompanyOutputCommunicator.empty()) return;

    while (true) {
        LOG("WORK")
        auto smartPolygon = m_PolygonsToSolve.pop();
        if (smartPolygon.m_Polygon == nullptr) {
            if (--m_CompanyCounter > 0) continue;
            unique_lock lock (g_SolverMutex);
            if (m_MinSolver.m_Solver && !m_MinSolver.m_Packs.empty()) {
                CEnhancedSolver solver (std::move(m_MinSolver));
                solver.solveAndMarkSolved();
            }
            if (m_CntSolver.m_Solver && !m_CntSolver.m_Packs.empty()) {
                CEnhancedSolver solver (std::move(m_CntSolver));
                solver.solveAndMarkSolved();
            }
            m_PolygonsToSolve.push(smartPolygon);
            break;
        }


        auto processProblem =
                [&](auto & polygon, CEnhancedSolver &solver, AProgtestSolver (*createSolver)()) {
            unique_lock lock (g_SolverMutex);
            if (!solver.m_Solver) solver = fetchSolver(createSolver);
            solver.m_Solver->addPolygon(polygon.m_Polygon);
            solver.m_Packs.push_back(polygon.m_EnhancedPack);
            if (!solver.m_Solver->hasFreeCapacity()) {
                CEnhancedSolver localSolver (std::move(solver));
                lock.unlock();
                localSolver.solveAndMarkSolved();
            }
        };
        smartPolygon.m_IsCnt ? processProblem(smartPolygon, m_CntSolver, createProgtestCntSolver)
                             : processProblem(smartPolygon, m_MinSolver, createProgtestMinSolver);
        SLEEP
    }

}

CEnhancedSolver COptimizer::fetchSolver(AProgtestSolver (*createSolver)()) {
    lock_guard lock(m_Mutex);
    CEnhancedSolver res;
    do {
        res = CEnhancedSolver(createSolver());
        SLEEP
    } while (!res.m_Solver->hasFreeCapacity());

    return res;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef __PROGTEST__

int main() {
    COptimizer optimizer;
    vector<ACompanyTest> companies;
    for (int i = 0; i < (10); ++i) {
        ACompanyTest company = std::make_shared<CCompanyTest>();
        optimizer.addCompany(company);
        companies.push_back(std::move(company));
    }

    optimizer.start(10);
    optimizer.stop();
    for (const auto &company: companies) {
        if (!company->allProcessed())
            throw std::logic_error("(some) problems were not correctly processsed");
    }
    return 0;
}

#endif /* __PROGTEST__ */
