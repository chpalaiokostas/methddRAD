export bam_to_sorted_bed
# turn bam files to a pseudo bed format
# the bed files are 1-based and not regular 0-based ones

"""
    bam_to_bed(filename::AbstractString)
creates a pseudo 1-based bed file 
"""    
function bam_to_bed_temp(bam::AbstractString)
    reader = open(BAM.Reader, bam)
    record = BAM.Record()
    output_name = string(split(basename(bam),".")[1], ".tmp")
    open(output_name, "w") do io
        while !eof(reader)
            empty!(record)
            read!(reader, record)
            if !BAM.ismapped(record) || !startswith(BAM.refname(record),"CM") || !(0 < BAM.templength(record) < 1000)
                continue
            end
            Chr=BAM.refname(record) 
            Start=BAM.position(record) # bed is 0-based but not changed here
            End=BAM.position(record) + BAM.templength(record)
            line = join([Chr,Start,End], "\t")
            println(io,line)
        end
    end
end

"""
    bam_to_sorted_bed(bam::AbstractString)
Converts a BAM file to a sorted, BED-like file using an external sort command.
"""
function bam_to_sorted_bed(bam::AbstractString)
    temp_name = string(split(basename(bam),".")[1], ".tmp")
    output_name = string(split(basename(bam),".")[1], ".bed")
    bam_to_bed_temp(bam)
    sort_command = `sort -k1,1 -k2,2n $temp_name`
    try
        run(pipeline(sort_command, stdout=output_name))
        println("✅ Successfully created sorted bed file: $output_name")
    catch e
        println("❌ Error during sorting: ", e)
    finally
        for file in readdir(".") 
            if endswith(file, ".tmp")
                rm(file)
            end    
        end
    end
end