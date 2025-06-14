export test_bootstrap, test_nobootstrap, compare_bootstrap
export test_reynaud2014

function test_bootstrap(model_sim::ParametricModel,
                        model_est::ParametricModel,
                        statistic::Type{<:Statistic},
                        A::Union{Interval, Area}
                        ; n_tests = 1000,
                        n_sims=1000)
   ps = zeros(n_tests)
   for i in 1:n_tests
       (i % 100 == 0) && println("Iteration $i")
       sim   = simulate(model_sim, A)
       ps[i] = fit_test(model_est, statistic, sim, n_sims=n_sims)
   end
   return ps
end

function test_nobootrap(model_sim::ParametricModel,
                        model_est::ParametricModel,
                        statistic::Type{<:Statistic},
                        A::Union{Interval, Area}
                        ; n_tests=1000,
                        n_sims=1000)
   ps = zeros(n_tests)
   for i in 1:n_tests
       (i % 100 == 0) && println("Iteration $i")
       sim   = simulate(model_sim, A)
       ps[i] = no_bootstrap(model_est, statistic, sim, n_sims=n_sims)
   end
   return ps
end

function compare_bootstrap(model_sim::ParametricModel,
                           model_est::ParametricModel,
                           statistic::Type{<:Statistic},
                           A::Union{Interval, Area},
                           ; n_tests=1000,
                           n_sims=1000)
    ps_boot  = zeros(n_tests)
    ps_nboot = zeros(n_tests)
    for i in 1:n_tests
        (i % 100 == 0) && println("Iteration $i")
        sim         = simulate(model_sim, A)
        ps_boot[i]  = fit_test(model_est, statistic, sim, n_sims=n_sims)
        ps_nboot[i] = no_bootstrap(model_est, statistic, sim, n_sims=n_sims)
    end
    return (ps_boot = ps_boot, ps_nboot = ps_nboot)
end

function test_reynaud2014(model_sim::ParametricModel,
                          model_est::ParametricModel,
                          A::Union{Interval, Area}
                          ; n_tests=1000)
   ps = zeros(n_tests)
   for i in 1:n_tests
       (i % 100 == 0) && println("Iteration $i")
       sim   = simulate(model_sim, A)
       ps[i] = reynaud2014(model_est, sim)
   end
   return ps
end
