
# Each itteration takes around 1 hour for 1million reads 32 threads

# Biggest improvement would be to rewrite RepeatExplorer in Rust and add GPU_Blast support
# Probably a 10x improvement in speed and efficiency
# 4x from GPU Blast alone and there is a lot of Blasting, 6x from writing in Rust vs Python and R

# Realistically right now, I can improve the script by
# Temporary File removal to prevent memory and file size bloat
# Check speed of Diamond vs blast
# Make one conda environment so I dont have to activate and deactive two seperate condas


# Would it be possible to first use the salamander repeat library to filter out the TE reads, and then use the remaining reads to build a new library. 
# This would potentially limit the specificity, but could reduce the complexity.