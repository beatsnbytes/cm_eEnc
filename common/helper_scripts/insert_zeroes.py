def insert_zeroes(input_file, output_file, k, l):
    with open(input_file, 'r') as infile:
        numbers = [int(line.strip()) for line in infile]

    result = []
    for i, num in enumerate(numbers):
        result.append(num)
        if (i + 1) % l == 0:
            result.extend([0] * k)

    with open(output_file, 'w') as outfile:
        for num in result:
            outfile.write(f"{num}\n")

if __name__ == "__main__":
    input_filename = "/home/vitis_wksp_212/decryption_cme/src/mceliece8192128/testbench/io_values/enc_4/pk_in.dat"  # Replace with the input filename
    output_filename = "/home/vitis_wksp_212/decryption_cme/src/mceliece8192128/testbench/io_values/enc_4/pk_in_alt_64.dat"  # Replace with the output filename
    k = 16  # Number of zeroes to insert
    l = 816  # Insert zeroes every l rows

    insert_zeroes(input_filename, output_filename, k, l)
