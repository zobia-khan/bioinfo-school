import random

def generate_random_read(length, gc_content=0.45, qual_profile="high"):
    """
    Generates a single random nucleotide read and corresponding quality scores.
    """
    bases = []
    quals = []
    
    # Base generation
    for _ in range(length):
        if random.random() < 0.02:  # 2% chance of ambiguous N
            bases.append("N")
        else:
            if random.random() < gc_content:
                bases.append(random.choice(["G", "C"]))
            else:
                bases.append(random.choice(["A", "T"]))
                
    # Quality score profile generation (Phred-33)
    # Phred score = ord(char) - 33
    # High quality: mostly Q30-Q40 (chr 63 to 73)
    # Low quality: declines towards the end, mostly Q10-Q20 (chr 43 to 53)
    for i in range(length):
        if qual_profile == "high":
            q = random.randint(30, 40)
        elif qual_profile == "declining":
            # Quality drops from Q38 at the start to Q12 at the end
            progress = i / length
            q = int(38 - progress * 26 + random.randint(-4, 4))
            q = max(0, min(41, q))
        else:  # low quality
            q = random.randint(8, 25)
            
        quals.append(chr(q + 33))
        
    return "".join(bases), "".join(quals)

def main():
    import os
    dirname = os.path.dirname("test.fastq")
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    
    # Generate 1000 mock reads
    random.seed(42)
    
    reads_to_generate = [
        # (count, length, gc_content, qual_profile, name_prefix)
        (400, 150, 0.45, "high", "READ_HQ"),
        (400, 150, 0.52, "declining", "READ_DEC"),
        (180, 100, 0.38, "low", "READ_LQ"),
        (20, 150, 0.75, "high", "READ_HIGH_GC")  # High GC spike
    ]
    
    with open("test.fastq", "w") as f:
        read_idx = 1
        for count, length, gc, profile, prefix in reads_to_generate:
            for _ in range(count):
                seq, qual = generate_random_read(length, gc, profile)
                
                # Introduce occasional Illumina TruSeq Adapter sequences at the end of some reads
                if prefix == "READ_DEC" and random.random() < 0.15:
                    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
                    # Pad adapter to be long enough
                    adapter_extended = adapter * 5
                    clip_pos = random.randint(80, 130)
                    seq = seq[:clip_pos] + adapter_extended[:length - clip_pos]
                    # Adjust quality scores for adapter to be medium
                    qual = qual[:clip_pos] + "".join(chr(random.randint(20, 30) + 33) for _ in range(length - clip_pos))
                
                f.write(f"@{prefix}_{read_idx} Description\n")
                f.write(f"{seq}\n")
                f.write("+\n")
                f.write(f"{qual}\n")
                read_idx += 1
                
    print(f"Generated {read_idx - 1} mock reads in 'test.fastq'.")

if __name__ == "__main__":
    main()
