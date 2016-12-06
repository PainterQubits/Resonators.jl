module Resonators
using InstrumentControl
using InstrumentControl.VNA
using InstrumentControl.E5071C

export search_sidecoupled

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
