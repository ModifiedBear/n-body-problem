# Importing
using CSV
using DataFrames
using Plots
using Random

function get_verlet(n::Int, m::Vector, G::Float64, t, dt::Float64, pos::Matrix, vel::Matrix)
    # initialize
    position = zeros(Float64, n,3,length(t))
    velocity = zeros(Float64, n,3,length(t))

    # initial pos and vel
    position[:,:,1] = pos
    velocity[:,:,1] = vel

    for ii=1:length(t) - 1
        # get dummy position
        posdummy = pos + vel * dt
        # get velocity
        vel = vel + 0.5 * get_acceleration(posdummy, m, G) * dt;
        # get position
        pos = posdummy + 0.5 * dt * vel;

        # append
        position[:,:,ii+1]= pos;
        velocity[:,:,ii+1]= vel;  
    end

    return (position, velocity)
end

function get_acceleration(position::Matrix, m::Vector, G::Float64)
    # Obtains the Newtonian acceleration for an object
    # -----------
    # Parameters:
    # -----------
    #     position: Nx3 array of positions
    #     M: Nx1 arrays of masses of said bodies
    #     G: Gravitational constant
    #     where N is  the number of bodies
    
    # adding a small error is required I don't know why
    E = 1E-10;  
    x = position[:,1,:];
    y = position[:,2,:];
    z = position[:,3,:];

    # following the equation from wikipedia, we have
    # r_j - r_i can be seen as a column vector - row 
    # vector (rows are "j", cols are "i" in a matrix 
    # M_ij, or rather M_ji since the operation is r_j
    # - r_i).
    
    # Get radius
    r = @.( ((x'- x)^2 + (y'- y)^2 + (z'- z)^2 + E^2)^0.5 );
    
    # Get cartesian accelerations
    # it HAS to be matrix multiplied with the masses
    ax = G * (x' .- x) ./ r.^3 * m;
    ay = G * (y' .- y) ./ r.^3 * m;
    az = G * (z' .- z) ./ r.^3 * m;
    
    # no need to specify Vector nor Type (Float64) as 
    # it already outputs it as such
    #A = Matrix{Float64}(undef, )
    A = [ax,ay,az];
    A = mapreduce(permutedims, vcat, A)

    return Matrix(transpose(A))

end

# Constants
# T0 = 0.
# TF = 600.
# DT = 0.1
# G = 1.
# 
# time = range(T0,TF,step=DT)
# 
# # read file
# data = CSV.read("bodies.csv", DataFrame)
# 
# # extract values
# pos0 = hcat(data.x, data.y, data.z)
# vel0 = hcat(data.xdot, data.ydot, data.zdot)
# 
# N = last(data.nbody)
# 
# M = data.m
# 
# 
# #pos, vel = get_verlet(N, M, G, time, DT, pos0, vel0);
# 
# position = zeros(Float64, N,3,length(time))
# velocity = zeros(Float64, N,3,length(time))
# 
# posvec, velvec = get_verlet(N, M, G, time, DT, pos0, vel0)

## Plotting

# base plot
# theme(:dark)
# gr(markersize=4)
# base = plot()
# 
# for ii = 1:N
#     
#     COLOR = RGB(rand(), 1,rand())
#     X = posvec[ii, 1, :]
#     Y = posvec[ii, 2, :]
#     Z = posvec[ii, 3, :]
#     
#     lab= "Body" * " " * string(ii)
#     
#     plot!(base, X, Y, Z,label=false, lc=COLOR)
#     scatter!(base, [last(X)], [last(Y)], [last(Z)], mc =COLOR, label=lab, msw=0)
#     scatter!(base, [first(X)], [first(Y)], [first(Z)], ma = 0., label=false)
# end
# current() # show current plot

function anim()

    # Constants
    T0 = 0.
    TF = 600.
    DT = 0.1
    G = 1.

    time = range(T0,TF,step=DT)

    # read file
    data = CSV.read("bodies.csv", DataFrame)

    # extract values
    pos0 = hcat(data.x, data.y, data.z)
    vel0 = hcat(data.xdot, data.ydot, data.zdot)

    N = last(data.nbody)

    M = data.m


    #pos, vel = get_verlet(N, M, G, time, DT, pos0, vel0);

    position = zeros(Float64, N,3,length(time))
    velocity = zeros(Float64, N,3,length(time))

    posvec, velvec = get_verlet(N, M, G, time, DT, pos0, vel0)
    

    COLOR = []
    for ii in 1:N
        COLOR = vcat(COLOR, RGB(rand(), 1,rand()))
    end

    plt = plot3d(1, dpi=300, legend = false)
    anim = @animate for ii in 1:100:length(time)
        # push!: do not clear prevous point
        #for jj ∈ 1:N

        jj = 1
        plt2 = plot!(plt, posvec[jj,1,1:ii], posvec[jj,2,1:ii], posvec[jj,3,1:ii], color=COLOR[jj],
            xlim=(-100,100), ylim=(-100, 100), zlim=(-100,100), alpha=max.((1:ii) .+ 800 .- ii, 0) / 800)
        #scatter!(plt, (posvec[jj,1,ii], posvec[jj,2,ii], posvec[jj,3,ii]), color=COLOR[jj])
        jj = 2
        plot!(plt, posvec[jj,1,1:ii], posvec[jj,2,1:ii], posvec[jj,3,1:ii], color=COLOR[jj],
            xlim=(-100,100), ylim=(-100, 100), zlim=(-100,100), alpha=max.((1:ii) .+ 800 .- ii, 0) / 800)
        #end
    end
    gif(anim, "body.gif", fps=60)
end
