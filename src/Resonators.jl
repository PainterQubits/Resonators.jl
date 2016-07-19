module Resonators

using PainterQB
using PainterQB.VNA
using PainterQB.E5071C
using DataFrames

import PainterQB.measure

export search_sidecoupled

function measure(x::VNA.FSweep)
    df = DataFrame()

    old_timeout = x.ins[Timeout]
    old_avgtrig = x.ins[AveragingTrigger]

    x.ins[Timeout] = 10000
    x.ins[NumTraces] = 3
    x.ins[VNA.Graphs] = [1 2; 1 3]
    x.ins[VNA.Format] = :PolarComplex
    x.ins[VNA.Format, 1, 2] = :LogMagnitude
    x.ins[VNA.Format, 1, 3] = :Phase
    x.ins[VNA.Parameter, 1, 2] = :S21
    x.ins[VNA.Parameter, 1, 3] = :S21
    x.ins[VNA.Parameter] = :S21
    x.ins[TriggerSource] = :BusTrigger

    if x.reject > 0
        old_avg = x.ins[Averaging]
        x.ins[Averaging] = false
        for i in 1:x.reject
            trig1(x.ins)
            sleep(sweeptime(x.ins))
            opc(x.ins)
        end
        x.ins[Averaging] = old_avg
    end

    x.ins[Averaging] && (x.ins[AveragingTrigger] = true)

    trig1(x.ins)
    sleep(sweeptime(x.ins))
    opc(x.ins)

    df[:f] = stimdata(x.ins)
    for s in [:S11, :S21]
        x.ins[VNA.Parameter] = s
        df[s] = VNA.data(x.ins, :PolarComplex)
    end

    autoscale(x.ins,1,1)
    autoscale(x.ins,1,2)
    autoscale(x.ins,1,3)

    # respect previous settings
    x.ins[Timeout] = old_timeout
    x.ins[AveragingTrigger] = old_avgtrig
    df
end

function search_sidecoupled(ins::InstrumentVNA, startf::Real, stopf::Real, cutoff=3, span=5e5)
    startf > stopf && error("Start frequency > stop frequency.")

    # Calculate frequency steps
    stepf = stopf/span
    numpts = length(startf:stepf:stopf)

    ins[TriggerSource] = :BusTrigger
    ins[FrequencyStart] = startf
    ins[FrequencyStop] = stopf
    ins[NumPoints] = numpts
    ins[Averaging] = false

    ins[VNA.Graphs] = [1]
    ins[NumTraces] = 1
    ins[VNA.Format, 1, 1] = :LogMagnitude
    ins[VNA.Parameter] = :S21
    ti = ins[SweepTime]

    trig1(ins)
    sleep(ti)
    opc(ins)

    # shotgun(ins, markers)
    markers = 1:nummarkers(ins)
    for m in markers
        ins[VNA.Marker, m] = true
        ins[VNA.MarkerX, m] = startf
    end

    peaks = Float64[]
    for (i,mi) in enumerate(markers)
        println(i)
        mx = search(ins, MarkerSearch(:RightPeak, 1, 1, mi, cutoff, VNA.Negative()))
        println(mx)
        isnan(mx) && break
        push!(peaks, mx)
        for mj in (i+1):nummarkers(ins)
            ins[VNA.MarkerX, mj] = mx+span
        end
    end

    for m in (length(peaks)+1):nummarkers(ins)
        ins[VNA.Marker, m] = false
    end

    # isnan(mx[1]) && (ins[VNA.Marker, 1] = false)
    # for m in 2:length(markers)
    #     (mx[m] == mx[m-1] || isnan(mx[m])) && (ins[VNA.Marker, m] = false)
    # end

    # fs = Float64[]
    # for m in markers
    #     if ins[VNA.Marker, m] == true
    #         count += 1
    #         push!(fs, mx[m])
    #     end
    # end

    # if count == 0
    if length(peaks) == 0
        warn("No peaks were found.")
    end

    autoscale(ins, 1, 1)
    peaks

    # count == 2 && configure(ins, Graphs, [1 2])
    # 3 <= count <= 4 && configure(ins, Graphs, [1 2; 3 4])
    # 5 <= count <= 6 && configure(ins, Graphs, [1 2 3; 4 5 6])
    # 7 <= count <= 8 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8])
    # count == 9 && configure(ins, Graphs, [1 2 3; 4 5 6; 7 8 9])
    # 10 <= count <= 12 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])
    # 13 <= count && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])
end



end # module
