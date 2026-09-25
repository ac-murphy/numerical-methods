#pragma once
#include <functional>

struct ode_state_1d
{
    float x, y;
};

struct ode_solution_1d
{
    std::vector<float> x_values, y_values;
};

struct wave_eq_initial_conditions_1d
{
    float c;
    std::pair<float, float> delta;
    std::pair<float, float> t_range;
    std::pair<float, float> x_range;
    std::pair<float, float> boundary_conditions;
    std::function<float(float)> initial_position;
    std::function<float(float)> initial_velocity;
};

struct wave_eq_solution_1d
{
    std::vector<float> t_values, x_values;
    std::vector<std::vector<float>> u_values;
};

inline ode_state_1d solve_RK4_step(const std::function<float(float, float)>& f, const float& xi, const float& yi, const float& dx)
{
    const float k1 = f(xi,             yi);
    const float k2 = f(xi + dx * 0.5f, yi + dx * k1 * 0.5f);
    const float k3 = f(xi + dx * 0.5f, yi + dx * k2 * 0.5f);
    const float k4 = f(xi + dx,        yi + dx * k3);

    return {
        xi + dx,
        yi + dx * (k1 + 2.0f * k2 + 2.0f * k3 + k4) / 6.0f,
    };
}
inline ode_solution_1d solve_RK4(const std::function<float(float, float)>& f, const float& x0, const float& y0, const float& x1, const float& dx)
{
    ode_solution_1d solution;
    solution.x_values.push_back(x0);
    solution.y_values.push_back(y0);

    while (solution.x_values.back() < x1)
    {
        const float xi = solution.x_values.back();
        const float yi = solution.y_values.back();

        const auto [x, y] = solve_RK4_step(f, xi, yi, dx);
        solution.x_values.push_back(x);
        solution.y_values.push_back(y);
    }

    return solution;
}

inline auto wave_eq_1d(const wave_eq_initial_conditions_1d& conditions)
{
    const auto& c = conditions.c;
    const auto& [dt, dx] = conditions.delta;
    const auto& [t_min, t_max] = conditions.t_range;
    const auto& [x_min, x_max] = conditions.x_range;
    const auto& [L, R] = conditions.boundary_conditions;
    const auto& f = conditions.initial_position;
    const auto& g = conditions.initial_velocity;

    // Calculate max index for time.
    size_t T = 0;
    for (float t = t_min; t < t_max; t += dt)
        ++T;

    // Calculate max index for x-axis.
    size_t X = 0;
    for (float x = x_min; x < x_max; x += dx)
        ++X;

    std::vector state(T, std::vector(X, 0.0f));

    // Set boundary conditions.
    for (size_t n = 0; n < T; ++n)
    {
        state[n][0] = L;
        state[n][X - 1] = R;
    }
    // Set initial conditions.
    for (size_t i = 0; i < X; ++i)
        state[0][i] = f(i * dx);
    for (size_t i = 1; i < X - 1; ++i)
        state[1][i] = f(i * dx) + ((c*c * dt*dt) / (2.0f * dx*dx)) * (state[0][i + 1] - 2.0f * state[0][i] + state[0][i - 1]);

    // Solve state.
    const float r = (c * dt) / dx;
    const float r2 = r * r;
    for (size_t n = 1; n < T - 1; ++n)
    {
        for (size_t i = 1; i < X - 1; ++i)
        {
            const float& x_prev = state[n - 1][i];
            const float& x_l = state[n][i - 1];
            const float& x_i = state[n][i];
            const float& x_r = state[n][i + 1];
            const float new_x_i = 2.0f * x_i - x_prev + r2 * (x_r - 2.0f * x_i + x_l);
            state[n + 1][i] = new_x_i;
        }
    }

    // Format solution.
    wave_eq_solution_1d solution;
    solution.u_values = state;
    for (size_t n = 0; n < T; ++n) solution.t_values.push_back(n * dt);
    for (size_t i = 0; i < X; ++i) solution.x_values.push_back(i * dx);

    // return solution;
    return solution;
}


