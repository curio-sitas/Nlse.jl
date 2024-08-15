using DifferentialEquations
const DE = DifferentialEquations

function laser!(du,u,p,s) 
    (ϵ, α, r, m, η) = p
    n = u[2]
    e =  u[1]
    du[1] = 0.5(1.0+α*1im)*n*e
    J = (1+m*cos(η*s))
    J = J*(J>=0)
    du[2] = ϵ*(r*J - (1+n)*(1+abs2(e)))

end

function integrate(plas, tseg, T, N)

(α, r, m, f) = plas

τphoton = 180e-12
τelec = 5e-12

ps = (τelec/τphoton,α, r ,m, f*τelec*2pi)


u0 = [0.1+0.0im, 0.01]
prob = ODEProblem(laser!, u0, (0.0, (tseg+T)/τelec), ps)
sol = DE.solve(prob, alg=RK4())

t = LinRange(tseg,tseg+T,N)

sol_n = sol(t./τelec)

return (sol_n[1,:],sol_n[2,:],t.-t[1])

end
