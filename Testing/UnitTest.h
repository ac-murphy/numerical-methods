#pragma once
#include <numbers>
#include "gtest/gtest.h"
#include "visualise.h"
#include "DifferentialEquations.h"

class UnitTest : public ::testing::Test
{
public:
    std::vector<float> sample_1d(const std::function<float(float)>& f, const std::vector<float>& x_values)
    {
        std::vector<float> y_values;
        y_values.reserve(x_values.size());
        for (const float& xi : x_values)
            y_values.push_back(f(xi));

        return y_values;
    }
};

#define ASSERT_DERIVATIVE_OF(df, f_expected, init_x, init_y, end_x, step) \
    const auto f_estimate = solve_RK4(df, init_x, init_y, end_x, step); \
    const auto y_expected = sample_1d(f_expected, f_estimate.x_values); \
    for (size_t i = 0; i < f_estimate.x_values.size(); ++i)      \
    {                                                            \
        ASSERT_NEAR(f_estimate.y_values[i], y_expected[i], 1e-5f); \
    }

TEST_F(UnitTest, RK4)
{
    auto df = [](const float& x, const float& y) { return -1.0f + 2.0f * x; };
    auto f = [](const float& x) { return -x + x * x; };
    ode_solution_1d solution = solve_RK4(df, 0.0f, f(0.0f), 1.0f, 0.01f);

    ASSERT_DERIVATIVE_OF(df, f, 0.0f, 0.0f, 1.0f, 0.01f);

    visualise vis;
    vis.init();
    vis.graph_1d(solution.x_values, solution.y_values);
    vis.graph_1d(solution.x_values, y_expected);
    vis.run();
}
TEST_F(UnitTest, FD_Wave_1D)
{
    const float L = 1.0f;
    const float h = 1.0f;
    wave_eq_initial_conditions_1d conditions;
    conditions.c = 5.0f;
    conditions.delta = { 0.001f, 0.01f };
    conditions.t_range = { 0.0f, 1.0f };
    conditions.x_range = { 0.0f, L };
    conditions.boundary_conditions = { 0.0f, 0.0f };
    conditions.initial_position = [&](const float& x){ return x > L * 0.5f ? (2.0f * h / L) * (L - x) : (2.0f * h / L) * x; };
    conditions.initial_velocity = [&](const float& x){ return 0.0f; };

    auto solution = wave_eq_1d(conditions);
    visualise vis;
    vis.init();
    vis.graph_1d_anim(solution.t_values, solution.x_values, solution.u_values);
    vis.run();
}
