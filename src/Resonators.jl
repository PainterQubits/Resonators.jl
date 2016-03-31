module Resonators

using PainterQB
using PainterQB.VNA
using PainterQB.E5071C
using DataFrames

import PainterQB.measure

export search_sidecoupled

function measure(x::VNA.FSweep)
    df = DataFrame()

    x.ins[NumTraces] = 3
    x.ins[VNA.Graphs] = [1 2; 1 3]
    x.ins[VNA.Format] = :PolarComplex
    x.ins[VNA.Format, 1, 2] = :LogMagnitude
    x.ins[VNA.Format, 1, 3] = :Phase
    x.ins[VNA.Parameter, 1, 2] = :S21
    x.ins[VNA.Parameter, 1, 3] = :S21
    x.ins[VNA.Parameter] = :S21
    ti = x.ins[SweepTime, 1]
    if (x.ins[Averaging] & x.ins[AveragingTrigger]) == true
        ti *= x.ins[AveragingFactor]
    end

    x.ins[TriggerSource] = :BusTrigger
    trig1(x.ins)
    sleep(ti+10)
    opc(x.ins)

    df[:f] = stimdata(x.ins)
    for s in [:S11, :S21]
        x.ins[VNA.Parameter] = s
        df[s] = VNA.data(x.ins, :PolarComplex)
    end

    autoscale(x.ins,1,1)
    autoscale(x.ins,1,2)
    autoscale(x.ins,1,3)

    df
end

function search_sidecoupled(ins::InstrumentVNA, startf::Real, stopf::Real, cutoff=3, span=5e5, markers=1:9)
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

    shotgun(ins, markers)
    mx = search(ins,
        [MarkerSearch(:RightPeak, 1, 1, m, cutoff, VNA.Negative()) for m in markers]...)

    isnan(mx[1]) && (ins[VNA.Marker, 1] = false)
    for m in 2:length(markers)
        (mx[m] == mx[m-1] || isnan(mx[m])) && (ins[VNA.Marker, m] = false)
    end

    count = 0
    fs = Float64[]
    for m in markers
        if ins[VNA.Marker, m] == true
            count += 1
            push!(fs, mx[m])
        end
    end

    if count == 0
        warning("No peaks were found.")
    end

    autoscale(ins, 1, 1)
    fs

    # count == 2 && configure(ins, Graphs, [1 2])
    # 3 <= count <= 4 && configure(ins, Graphs, [1 2; 3 4])
    # 5 <= count <= 6 && configure(ins, Graphs, [1 2 3; 4 5 6])
    # 7 <= count <= 8 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8])
    # count == 9 && configure(ins, Graphs, [1 2 3; 4 5 6; 7 8 9])
    # 10 <= count <= 12 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])
    # 13 <= count && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])
end



end # module
