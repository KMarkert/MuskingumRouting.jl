module MuskingumRouting

function rout(inflow, k, x, dt)
    C0,C1,C2 = calc_cs(k,x,dt)

    len = length(inflow)
    outflow = zeros(Float64,len)
    outflow[1] = inflow[1]

    for i in 1:(len-1)
        ni = inflow[i+1]
        ti = inflow[i]
        to = outflow[i]

        outflow[i+1] = (C0 * ni) + (C1 * ti) + (C2 * to)
    end

    return outflow
end

function calc_cs(k, x, dt)
    Δt = (0.5 * dt)
    invx = k * (1 - x)
    C0 = ((-k * x) + Δt) / (invx + Δt)
    C1 = ((k * x) + Δt) / (invx + Δt)
    C2 = (invx - Δt) / (invx + Δt)

    return C0, C1, C2
end

function calc_params(inflow,outflow,dt,init_storage)

    s = calc_storage(inflow,outflow,dt,init_storage)

    o2 = sum(outflow .^ 2)
    si = sum(s*inflow)
end

function calc_storage(inflow,outflow,dt,init_storage)
    len = length(inflow)
    storage = zeros(Float64,len)
    storage[1] = init_storage

    # Δt = (0.5 * dt)

    for i in 1:(len-1)
        storage[i+1] = storage[i] + (dt * (inflow[i+1] + inflow[i]))/2 - (dt* (outflow[i+1] + outflow[i]))/2
    end

    return storage
end

end # module
