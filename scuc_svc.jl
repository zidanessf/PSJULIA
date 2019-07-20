using JuMP,PowerModels,DataFrames,CSV,CPLEX,Random,Suppressor
using Redis,DataStructures,JSON,Base64
function process(conn::RedisConnection)
    while true
        try
            if !isempty(qsub)
                REQ = dequeue!(qsub)
                d = JSON.parse(REQ)
                confidence = d["confidence"]
                k = d["windCur_weight"]/5
                @info("收到计算请求，置信度为：$confidence\n风电消纳指标权重为:$k")
                println("正在解析算例文件...")
                caseBytes = base64decode(d["case"])
                filename = d["ID"]*".m"
                io = open(filename,"w")
                write(io,caseBytes)
                close(io)
                println("算例文件储存为:"*filename*"")
                @info("开始计算经济运行域:\n")
                resp = GetEconomicalRegion(d,filename)
                publish(conn,"opt_response",resp)
            end
        catch exception
            println(exception)
            if isa(exception,InterruptException)
                break
            end
        end
        sleep(0.5)
    end
end
function GetEconomicalRegion(REQ,casefile)
    RESP = Dict("objective"=>0)
    k = REQ["windCur_weight"]/5
    m = JuMP.Model(with_optimizer(CPLEX.Optimizer))
    silence()
    case = PowerModels.parse_file(casefile)
    pm = build_model(case, ACPPowerModel,
    PowerModels.post_opf)
    ref = pm.ref[:nw][0]
    for gen in keys(ref[:gen])
        if ref[:gen][gen]["pmin"] == 0
            ref[:gen][gen]["pmin"] = 0.1 * ref[:gen][gen]["pmax"]
        end
        if ref[:gen][gen]["pmax"] == 0
            delete!(ref[:gen],gen)
        end
    end
    T = 96
    aggr = []
    gens = setdiff([gen for gen in keys(ref[:gen])],aggr)
    pmax = maximum([ref[:gen][gen]["pmax"] for gen in keys(ref[:gen])])
    hydro = Dict("pmax_pump"=>pmax/6,"pmin_pump"=>0.1,"pmax_gen"=>pmax/6,"pmin_gen"=>0.1,"e2w"=>2.5,"w2e"=>1.6,"C_max"=>154.900,"C_min"=>37.820,"C_int"=>97.820)
    ref[:hydro] = Dict([(i,hydro) for i in range(1,6)])
    hydros = [gen for gen in keys(ref[:hydro])]
    price = @variable(m,[1:T],lower_bound=0.5,upper_bound=1.5)
    pw = @variable(m,[x in gens,1:T])
    hydro_gen = @variable(m,[x in hydros,1:T])
    hydro_pump = @variable(m,[x in hydros,1:T])
    hydro_v = @variable(m,[x in hydros,1:T])
    hydro_W = @variable(m,[x in hydros,1:T])
    u_gen = @variable(m,[x in hydros,1:T])
    u_pump = @variable(m,[x in hydros,1:T])
    pv = @variable(m,[x in gens],lower_bound=0)
    windv = @variable(m,[1:T],lower_bound=0)
    pub = @variable(m,[x in gens,1:T])
    plb = @variable(m,[x in gens,1:T])
    st = @variable(m,[x in gens,1:T],Bin)
    st_z =@variable(m,[x in gens,1:T],Bin)
    st_y =@variable(m,[x in gens,1:T],Bin)
    load = sum([ref[:load][x]["pd"] for x in keys(ref[:load])])
    Random.seed!(1234)
    # load_real = load + 0.2*rand(Float64,T)*load
    load_data = REQ["system_load"]
    load_real = load_data * load/maximum(load_data)
    wind_data = REQ["WIND"][1]##TODO
    ramp = 0.5
    r = 0.3 * load/maximum(wind_data["upper_bound"])
    windmax = r * wind_data["upper_bound"]
    windmid = r * (wind_data["upper_bound"] + wind_data["lower_bound"])/2
    windmin = r * wind_data["lower_bound"]
    pmax = [ref[:gen][x]["pmax"] for x in aggr]
    for t in 1:T
        for gen in hydros
            @constraint(m,hydro_gen[gen,t] >= u_gen[gen,t] * ref[:hydro][gen]["pmin_gen"])
            @constraint(m,hydro_gen[gen,t] <= u_gen[gen,t] * ref[:hydro][gen]["pmax_gen"])
            @constraint(m,hydro_pump[gen,t] >= u_pump[gen,t] * ref[:hydro][gen]["pmin_pump"])
            @constraint(m,hydro_pump[gen,t] <= u_pump[gen,t] * ref[:hydro][gen]["pmax_pump"])
            @constraint(m,u_pump[gen,t] + u_gen[gen,t] <= 1)
            @constraint(m,hydro_gen[gen,t] - hydro_pump[gen,t] - (windmax[t] - windmid[t])*hydro_v[gen,t] >=  - u_pump[gen,t] * ref[:hydro][gen]["pmax_pump"])
            @constraint(m,hydro_gen[gen,t] - hydro_pump[gen,t] - (windmax[t] - windmid[t])*hydro_v[gen,t] <=  u_gen[gen,t] * ref[:hydro][gen]["pmax_gen"])
            if t >=2
                @constraint(m,(windmax[t] - windmid[t])*hydro_v[gen,t] - (windmin[t-1] - windmid[t-1])*hydro_v[gen,t] <= 0.2 * ref[:hydro][gen]["pmax_gen"])
                @constraint(m,(windmax[t] - windmid[t])*hydro_v[gen,t] - (windmax[t-1] - windmid[t-1])*hydro_v[gen,t] >= -0.2 * ref[:hydro][gen]["pmax_gen"])
            end
            if t == 1
                @constraint(m,hydro_W[gen,t] == ref[:hydro][gen]["C_int"])
            else
                @constraint(m,hydro_W[gen,t] == hydro_W[gen,t-1] + ref[:hydro][gen]["e2w"]*hydro_pump[gen,t] - ref[:hydro][gen]["e2w"]*hydro_gen[gen,t])
            end
            @constraint(m,hydro_W[gen,t]<=ref[:hydro][gen]["C_max"])
            @constraint(m,hydro_W[gen,t]>=ref[:hydro][gen]["C_min"])
            if t == 24
                @constraint(m,hydro_W[gen,t]==ref[:hydro][gen]["C_int"])
            end
        end
        for gen in gens
            @constraint(m,pub[gen,t] <= st[gen,t]*ref[:gen][gen]["pmax"])
            @constraint(m,plb[gen,t] >= st[gen,t]*ref[:gen][gen]["pmin"])
            @constraint(m,pw[gen,t] - (windmin[t] - windmid[t])*pv[gen]  <= pub[gen,t])
            @constraint(m,pw[gen,t] - (windmax[t] - windmid[t])*pv[gen] >= plb[gen,t])
            if t >=2
                @constraint(m,st[gen,t] == st[gen,t-1] + st_y[gen,t] - st_z[gen,t])
                @constraint(m,st_y[gen,t] + st_z[gen,t] <= 1)
                @constraint(m,pub[gen,t] - plb[gen,t-1] <= ramp*ref[:gen][gen]["pmax"]) # add adjustable var
                @constraint(m,plb[gen,t-1] - pub[gen,t] <= ramp*ref[:gen][gen]["pmax"]) # add adjustable var
                @constraint(m,sum(st_y[gen,k] for k in max(t-4,1):t)<=st[gen,t])
                @constraint(m,sum(st_z[gen,k] for k in max(t-4,1):t)<=1 - st[gen,t])
            else
                @constraint(m,st[gen,t] == st_y[gen,t] - st_z[gen,t])
                @constraint(m,st_y[gen,t] + st_z[gen,t] <= 1)
            end
        end
        for x in aggr
            @constraint(m, ref[:gen][x]["pmax"] - (windmin[t] - windmid[t])*aggr_v[x]  <= 1.2 * ref[:gen][x]["pmax"])
            @constraint(m, ref[:gen][x]["pmax"] - (windmax[t] - windmid[t])*aggr_v[x] >= 0.5 * ref[:gen][x]["pmax"])
        end
        # 不可调功率平衡
        if length(aggr) > 0
            @constraint(m,sum(pw[gen,t] for gen in gens)
            - (2 - price[t]) * sum([ref[:gen][x]["pmax"] for x in aggr])
            + windmid[t]
            + sum(hydro_gen[gen,t] for gen in hydros) == load_real[t] + sum(hydro_pump[gen,t] for gen in hydros))
        else
            @constraint(m,sum(pw[gen,t] for gen in gens)
            + windmid[t]
            + sum(hydro_gen[gen,t] for gen in hydros)== load_real[t] + sum(hydro_pump[gen,t] for gen in hydros))
        end
        # 可调功率平衡
        @constraint(m,sum(pv[gen] for gen in gens) + sum(hydro_v[gen,t] for gen in hydros) + sum(aggr_v[gen] for gen in aggr) + windv[t] == 1)
    end
    startup_cost = sum(ref[:gen][gen]["startup"]*st_y[gen,t] for gen in gens,t in 1:T-1)
    fuel_cost = sum(ref[:gen][gen]["cost"][2]*pw[gen,t] for gen in gens,t in 1:T) #+ ref[:gen][gen]["cost"][1]*pw[gen,t]^2
    run_cost = 0
    try
        run_cost = sum(ref[:gen][gen]["cost"][3] * st[gen,t] for gen in gens,t in 1:T) # May slow down the program
    catch exception
    end
    agc_cost = sum(0.1*hydro_v[gen,t]^2*(windmax[t]-windmid[t])^2 for gen in hydros,t in 1:T)
    interval = 100*sum(pub[gen,t] - plb[gen,t] for gen in gens,t in 1:T)
    CO2Emission = sum(900+(100000/ref[:gen][gen]["cost"][1])*pw[gen,t] for gen in gens,t in 1:T)
    windCur = 100*sum(windv[t]*(windmax[t]-windmin[t]) for t in 1:T)
    @objective(m,Min,startup_cost+ run_cost+ interval + fuel_cost + REQ["windCur_weight"]*windCur)
    try
        optimize!(m)
    catch e
        println(e)
    end
    result = Dict()
    result["power"] = value.(pw)
    result["on_off"] = value.(st)
    result["power_ub"] = value.(pub)
    result["power_lb"] = value.(plb)
    result["factor"] = value.(pv)
    nonagc = DataFrame()
    agc = DataFrame()
    wind = DataFrame()
    for gen in gens
        if value.(pv[gen]) == 0
            nonagc[Symbol("gen",gen," state")] = value.(st[gen,:])
            # nonagc[Symbol("gen",gen," on")] = value.(st_y[gen,:])
            # nonagc[Symbol("gen",gen," off")] = value.(st_z[gen,:])
            nonagc[Symbol("gen",gen," power")] = value.(pw[gen,:])
        else
            agc[Symbol("gen",gen," state")] = value.(st[gen,:])
            agc[Symbol("gen",gen," power")] = value.(pw[gen,:])
            agc[Symbol("gen",gen," upperbound")] = value.(pub[gen,:])
            agc[Symbol("gen",gen," lowerbound")] = value.(plb[gen,:])
        end
    end
    for gen in hydros
        agc[Symbol("hydro",gen,"power")] = Array(value.(hydro_gen[gen,:])) -  Array(value.(hydro_pump[gen,:]))
        agc[Symbol("hydro",gen," upperbound")] =  Array(value.(hydro_gen[gen,:])) -   Array(value.(hydro_pump[gen,:])) - (windmin - windmid).*[value(hydro_v[gen,t]) for t in 1:T]
        agc[Symbol("hydro",gen," lowerbound")] =  Array(value.(hydro_gen[gen,:])) -  Array(value.(hydro_pump[gen,:])) - (windmax - windmid).*[value(hydro_v[gen,t]) for t in 1:T]
    end
    wind[:min] = windmin
    wind[:mean] = windmid
    wind[:max] = windmax
    wind[:price] = value.(price)
    #publish response
    resp = Dict("ID"=>REQ["ID"],"agc"=>"","non_agc"=>"","wind"=>"")
    io = IOBuffer()
    CSV.write(io,agc)
    resp["agc"] = base64encode(String(take!(io)))
    CSV.write(io,nonagc)
    resp["nonagc"] = base64encode(String(take!(io)))
    CSV.write(io,wind)
    resp["wind"] = base64encode(String(take!(io)))
    resp_json = JSON.json(resp)
    #end publish
    cost = value.(startup_cost + fuel_cost + run_cost)
    CO2 = value.(CO2Emission)
    windPercent = sum(1 - value.(windv[t]) for t in 1:T)/T
    io = open("logfile.txt","a")
    write(io,string("k = ","$k","\n"))
    write(io,string("运行成本：","$cost","\n"))
    write(io,string("风电消纳率：","$windPercent","\n\n"))
    close(io)
    resp_json
end
qsub = Queue{Any}()
f(q::Queue{Any},y) = begin
    enqueue!(q,y)
end
conn = RedisConnection()
sub = open_subscription(conn)
subscribe(sub,"opt_requests",y->f(qsub,y))
@info("connection established")
# implementations of our test servies
t = @async process(conn)
@info("listening...")
#Base.throwto(t, InterruptException())
