using JuMP,PowerModels,DataFrames,Gurobi,CSV,Random,Suppressor
using Redis,DataStructures,JSON,Base64
const UT,DT,ramp = 16,16,0.25*0.25
function process(conn::RedisConnection)
    try
        while true
            if !isempty(qsub)
                REQ = dequeue!(qsub)
                d = JSON.parse(REQ)
                confidence = d["confidence"]
                k = d["windCur_weight"]
                @info("收到计算请求，置信度为：$confidence\n风电消纳指标权重为:$k")
                println("正在解析算例文件...")
                caseBytes = d["case"]
                filename = d["ID"]*".m"
                io = open(filename,"w")
                @info("file opened")
                write(io,caseBytes)
                @info("file written")
                close(io)
                @info("算例文件储存为:"*filename*"")
                @info("开始计算经济运行域:\n")
                pre_calc = makeWindScheme()
                resp = getEconomicalRegion(d,filename)
                publish(conn,"opt_response",resp)
            end
            sleep(0.5)
        end
    catch exception
        @error(exception)
        @error("出现异常，程序关断中...")
        @info("disconnecting Redis")
        unsubscribe(sub,"opt_requests")
        disconnect(conn)
        @info("程序结束")
        throw(exception)
    end
end
function makeWindScheme(ref::Dict,load,wind)

function get_executvie_summary(ref::Dict,load::Array,wind::Array)
    Dict("ngen" => length(keys(ref[:gen])),
    "nbus" => length(keys(ref[:bus])),
    "nbranch" => length(keys(ref[:branch])),
    "totalPmax" => sum(ref[:gen][gen]["pmax"] for gen in keys(ref[:gen])),
    "totalLoadMax" => maximum(load),
    "totalLoadMin" => minimum(load),
    "windPercentMax" => 0.3,
    "totalWindMax" => maximum(wind),
    "totalWindMin" => minimum(wind))
end
function getEconomicalRegion(REQ,casefile)
    RESP = Dict("objective"=>0)
    if REQ["windCur_weight"] < 0
        @warn("风电消纳权重小于0，重置为0")
    elseif REQ["windCur_weight"] <=100
        k_wind = 3 * REQ["windCur_weight"]
    else
        k_wind = 3
        @warn("风电消纳权重大于1，重置为1")
    end
    m = JuMP.direct_model(Gurobi.Optimizer(Presolve=1,Crossover=0,MIPGap=1e-3))
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
    gens_param_lost = [gen for gen in keys(ref[:gen]) if length(ref[:gen][gen]["cost"])<2]
    @warn("机组$gens_param_lost 参数缺失")
    gens = setdiff([gen for gen in keys(ref[:gen])],aggr)
    gens = setdiff(gens,gens_param_lost)
    pmax_dict = [(gen,ref[:gen][gen]["pmax"]) for gen in gens]
    gen_top_sort_by_pmax = Base.sort(pmax_dict,by=x->x[2],rev=true)[1:min(convert(Int,floor((length(gens)/2))),30)]
    gen_top_sort_by_pmax = [x[1] for x in gen_top_sort_by_pmax]
    pmax = maximum([ref[:gen][gen]["pmax"] for gen in keys(ref[:gen])])
    hydro = Dict("pmax_pump"=>pmax/6,"pmin_pump"=>0.1,"pmax_gen"=>pmax/6,"pmin_gen"=>0.1,"e2w"=>2.5,"w2e"=>1.6,"C_max"=>154.900,"C_min"=>37.820,"C_int"=>97.820,"om_cost"=>0.05)
    ref[:hydro] = Dict([(i,hydro) for i in range(1,6)])
    hydros = [gen for gen in keys(ref[:hydro])]
    price = @variable(m,[1:T],lower_bound=0.5,upper_bound=1.5)
    pw = @variable(m,[x in gens,1:T],lower_bound=0)
    hydro_gen = @variable(m,[x in hydros,1:T],lower_bound=0)
    hydro_pump = @variable(m,[x in hydros,1:T],lower_bound=0)
    hydro_v = @variable(m,[x in hydros,1:T],lower_bound=0,upper_bound=1)
    hydro_W = @variable(m,[x in hydros,1:T])
    # u_gen = @variable(m,[x in hydros,1:T],Bin)
    # u_pump = @variable(m,[x in hydros,1:T],Bin)
    pv = @variable(m,[x in gens],lower_bound=0,upper_bound=1)
    windv = @variable(m,[1:T],lower_bound=0,upper_bound=1)
    pub = @variable(m,[x in gens,1:T],lower_bound=0)
    plb = @variable(m,[x in gens,1:T],lower_bound=0)
    st = @variable(m,[x in gens,1:T],Bin)
    st_aux = @variable(m,[x in gens, 1:T],lower_bound=0)
    # st_z =@variable(m,[x in gens,1:T],Bin)
    # st_y =@variable(m,[x in gens,1:T],Bin)
    # Set Bounds
    for t in 1:T
        for gen in gens
            set_lower_bound(pw[gen,t],0)
            set_upper_bound(pw[gen,t],ref[:gen][gen]["pmax"])
            set_upper_bound(pub[gen,t],ref[:gen][gen]["pmax"])
            set_upper_bound(plb[gen,t],ref[:gen][gen]["pmax"])
        end
        for gen in hydros
            set_upper_bound(hydro_gen[gen,t],ref[:hydro][gen]["pmax_gen"])
            set_upper_bound(hydro_pump[gen,t],ref[:hydro][gen]["pmax_pump"])
            set_lower_bound(hydro_W[gen,t],ref[:hydro][gen]["C_min"])
            set_upper_bound(hydro_W[gen,t],ref[:hydro][gen]["C_max"])
        end
    end
    for gen in setdiff(gens,gen_top_sort_by_pmax)
        set_upper_bound(pv[gen],0)
    end
    for gen in gen_top_sort_by_pmax
        set_lower_bound(pv[gen],0.5*1/length(gen_top_sort_by_pmax))
    end
    # load = sum([ref[:load][x]["pd"] for x in keys(ref[:load])])
    load = 0.8*sum(ref[:gen][gen]["pmax"] for gen in gens)
    Random.seed!(1234)
    # load_real = load + 0.2*rand(Float64,T)*load
    load_data = REQ["system_load"]
    load_real = load_data * load/maximum(load_data)
    wind_data = REQ["WIND"][1]##TODO
    r = 0.4 * load/maximum(wind_data["upper_bound"])
    windmax = r * wind_data["upper_bound"]
    windmid = r * (wind_data["upper_bound"] + wind_data["lower_bound"])/2
    windmin = r * wind_data["lower_bound"]
    pmax = [ref[:gen][x]["pmax"] for x in aggr]
    case_info::Dict = get_executvie_summary(ref,load_real,windmid)
    for t in 1:T
        for gen in hydros
            # @constraint(m,hydro_gen[gen,t] >= u_gen[gen,t] * ref[:hydro][gen]["pmin_gen"])
            # @constraint(m,hydro_gen[gen,t] <= u_gen[gen,t] * ref[:hydro][gen]["pmax_gen"])
            # @constraint(m,hydro_pump[gen,t] >= u_pump[gen,t] * ref[:hydro][gen]["pmin_pump"])
            # @constraint(m,hydro_pump[gen,t] <= u_pump[gen,t] * ref[:hydro][gen]["pmax_pump"])
            # @constraint(m,u_pump[gen,t] + u_gen[gen,t] <= 1)
            @constraint(m,hydro_gen[gen,t] - hydro_pump[gen,t] - (windmax[t] - windmid[t])*hydro_v[gen,t] >=  - ref[:hydro][gen]["pmax_pump"])
            @constraint(m,hydro_gen[gen,t] - hydro_pump[gen,t] - (windmax[t] - windmid[t])*hydro_v[gen,t] <=  ref[:hydro][gen]["pmax_gen"])
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
            #@constraint(m,plb[gen,t] >= st[gen,t]*ref[:gen][gen]["pmin"])
            @constraint(m,pw[gen,t] - (windmin[t] - windmid[t])*pv[gen]  <= pub[gen,t])
            @constraint(m,pw[gen,t] - (windmax[t] - windmid[t])*pv[gen] >= plb[gen,t])
            if t >=2
                # @constraint(m,st[gen,t] == st[gen,t-1] + st_y[gen,t] - st_z[gen,t])
                # @constraint(m,st_y[gen,t] + st_z[gen,t] <= 1)
                @constraint(m,pub[gen,t] - plb[gen,t-1] <= ramp*ref[:gen][gen]["pmax"]) # add adjustable var
                @constraint(m,plb[gen,t-1] - pub[gen,t] <= ramp*ref[:gen][gen]["pmax"]) # add adjustable var
                # @constraint(m,sum(st_y[gen,k] for k in max(t-4,1):t)<=st[gen,t])
                # @constraint(m,sum(st_z[gen,k] for k in max(t-4,1):t)<=1 - st[gen,t])
                if t <= T - UT + 1
                    @constraint(m,sum(st[gen,k] for k in t:t+UT-1) >= UT*(st[gen,t] - st[gen,t-1]))
                else
                    @constraint(m,sum(st[gen,k] - (st[gen,t] - st[gen,t-1]) for k in t:T) >= 0)
                end
                if t <= T - DT + 1
                    @constraint(m,sum(1 - st[gen,k] for k in t:t+DT-1) >= DT*(st[gen,t-1] - st[gen,t]))
                else
                    @constraint(m,sum(1 - st[gen,k] - (st[gen,t-1] - st[gen,t]) for k in t:T) >= 0)
                end
                @constraint(m,st_aux[gen,t] >= st[gen,t] - st[gen,t-1])
            else
                #TODO add some initial conditions
                #@constraint(m,st[gen,t] == st_y[gen,t] - st_z[gen,t])
                # @constraint(m,st_y[gen,t] + st_z[gen,t] <= 0)
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
            + windmid[t] - wind_Cur[t]
            + sum(hydro_gen[gen,t] for gen in hydros) == load_real[t] + sum(hydro_pump[gen,t] for gen in hydros))
        else
            @constraint(m,sum(pw[gen,t] for gen in gens)
            + windmid[t]
            + sum(hydro_gen[gen,t] for gen in hydros)== load_real[t] + sum(hydro_pump[gen,t] for gen in hydros))
        end
        # 可调功率平衡
        @constraint(m,sum(pv[gen] for gen in gens) + sum(hydro_v[gen,t] for gen in hydros) + sum(aggr_v[gen] for gen in aggr) + windv[t] == 1)
        @info("constraints of time$t constructed")
    end
    #@constraint(m,lb_to_relax[gen=gens,t=1:T],plb[gen,t] >= st[gen,t]*ref[:gen][gen]["pmin"])
    startup_cost = @expression(m,sum(ref[:gen][gen]["startup"]*st_aux[gen,t] for gen in gens,t in 1:T))
    @info("startup_cost constructed")
    fuel_cost = @expression(m,sum((ref[:gen][gen]["cost"][end-1] + index/length(gens))*(pw[gen,t]) for (index,gen) in enumerate(gens),t in 1:T)) #+ ref[:gen][gen]["cost"][1]*pw[gen,t]^2
    @info("fuel_cost constructed")
    hydro_cost = @expression(m,sum(ref[:hydro][gen]["om_cost"]*(hydro_pump[gen,t] + hydro_gen[gen,t]) for gen in hydros,t in 1:T))
    @info("hydro_cost constructed")
    run_cost = 0
    try
        run_cost = @expression(m,sum(ref[:gen][gen]["cost"][3] * st[gen,t] for gen in gens,t in 1:T))# May slow down the program
    catch exception
    end
    @info("runcost constructed")
    agc_cost = @expression(m,sum(0.1*hydro_v[gen,t]^2*(windmax[t]-windmid[t])^2 for gen in hydros,t in 1:T))
    @info("agc_cost constructed")
    interval = @expression(m,10*sum(pub[gen,t] - plb[gen,t] for gen in gens,t in 1:T))
    @info("interval_cost constructed")
    CO2Emission = @expression(m,sum(900+(100000/ref[:gen][gen]["cost"][1])*pw[gen,t] for gen in gens,t in 1:T))
    @info("CO2Emission constructed")
    windCur = @expression(m,100*sum(windv[t]*(windmax[t]-windmin[t]) for t in 1:T))
    @info("windCur constructed")
    @objective(m,Min,startup_cost+ run_cost+ interval + fuel_cost + hydro_cost+ k_wind*windCur)
    @info("objective constructed")
    @info("建模完毕，首先进行松弛计算...")
    # 机组开停机状态变量固定为1
    try
        for t in 1:T
            for gen in gens
                unset_binary(st[gen,t])
                # unset_binary(st_y[gen,t])
                # unset_binary(st_z[gen,t])
                fix(st[gen,t],1)
                # fix(st_y[gen,t],0)
                # fix(st_z[gen,t],0)
            end
        end
        optimize!(m)
    catch e
        println(e)
    end
    @info("松弛计算完毕，缩小机组启停寻优范围...")
    num_have_to_be_on = 0
    num_have_to_be_off = 0
    for gen in gens
        on_flag = true
        off_flag = true
        for t in 1:T
            if value(pw[gen,t]) < 0.65 * ref[:gen][gen]["pmax"]
                on_flag = false
                break
            end
        end
        for t in 1:T
            if value(pw[gen,t]) > 0.01 * ref[:gen][gen]["pmax"]
                off_flag = false
                break
            end
        end
        if on_flag
            nothing
            num_have_to_be_on += 1
        # elseif off_flag
        #     num_have_to_be_off += 1
        #     for t in 1:T
        #         fix(st[gen,t],0)
        #     end
        else
            for t in 1:T
                unfix(st[gen,t])
                # unfix(st_z[gen,t])
                # unfix(st_y[gen,t])
                set_binary(st[gen,t])
                # set_binary(st_y[gen,t])
                # set_binary(st_z[gen,t])
            end
        end
    end
    @info("必开机组数量为$num_have_to_be_on")
    @info("必关机组数量为$num_have_to_be_off")
    @info("开始完整计算...")
    @constraint(m,lb_to_relax[gen=gens,t=1:T],plb[gen,t] >= st[gen,t]*ref[:gen][gen]["pmin"])
    optimize!(m)
    result = Dict()
    result["power"] = value.(pw)
    result["on_off"] = value.(st)
    result["power_ub"] = value.(pub)
    result["power_lb"] = value.(plb)
    result["factor"] = value.(pv)
    nonagc = DataFrame()
    agc = DataFrame()
    wind = DataFrame()
    sort_aux = [(x,(sum(value.(pub[x,:])) - sum(value.(plb[x,:])))/(ref[:gen][x]["pmax"]+0.001)) for x in gens]
    gens_sorted_by_interval = Base.sort(sort_aux,by=x->x[2],rev=true)
    gens_sorted_by_interval = [x[1] for x in gens_sorted_by_interval]
    for gen in gens_sorted_by_interval
        if value.(pv[gen]) == 0
            if sum(value.(pw[gen,:])) <= 0.01
                nonagc[Symbol("ALWAYS_OFF_gen",gen," state")] = value.(st[gen,:])
                nonagc[Symbol("ALWAYS_OFF_gen",gen," power")] = value.(pw[gen,:])
            else
                nonagc[Symbol("gen",gen," state")] = value.(st[gen,:])
                nonagc[Symbol("gen",gen," power")] = value.(pw[gen,:])
            end
        else
            agc[Symbol("gen",gen," state")] = value.(st[gen,:])
            agc[Symbol("gen",gen," power")] = value.(pw[gen,:])
            agc[Symbol("gen",gen," upperbound")] = value.(pub[gen,:])
            agc[Symbol("gen",gen," lowerbound")] = value.(plb[gen,:])
        end
    end
    for gen in hydros
        agc[Symbol("hydro",gen," power")] = Array(value.(hydro_gen[gen,:])) -  Array(value.(hydro_pump[gen,:]))
        agc[Symbol("hydro",gen," upperbound")] =  Array(value.(hydro_gen[gen,:])) -   Array(value.(hydro_pump[gen,:])) - (windmin - windmid).*[value(hydro_v[gen,t]) for t in 1:T]
        agc[Symbol("hydro",gen," lowerbound")] =  Array(value.(hydro_gen[gen,:])) -  Array(value.(hydro_pump[gen,:])) - (windmax - windmid).*[value(hydro_v[gen,t]) for t in 1:T]
    end
    wind[:min] = windmin
    wind[:mean] = windmid
    wind[:max] = windmax
    wind[:price] = value.(price)
    #publish response
    cost = value.(startup_cost + fuel_cost + run_cost)
    CO2 = value.(CO2Emission)
    windPercent = sum(1 - value.(windv[t]) for t in 1:T)/T
    resp = Dict("ID"=>REQ["ID"],"agc"=>"","non_agc"=>"","wind"=>"",
    "k_wind"=>k_wind,"cost"=>cost,"wind_acception"=>windPercent,
    "case_info"=>case_info)
    io = IOBuffer()
    CSV.write(io,agc)
    resp["agc"] = base64encode(String(take!(io)))
    CSV.write(io,nonagc)
    resp["nonagc"] = base64encode(String(take!(io)))
    CSV.write(io,wind)
    resp["wind"] = base64encode(String(take!(io)))
    resp_json = JSON.json(resp)
    #end publish
    io = open("logfile.txt","a")
    write(io,string("k = ","$k_wind","\n"))
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

while true
   sleep(1)
end
#Base.throwto(t, InterruptException())
