function mat2bin(s,bin)
[nt nh]=size(s);
fd=fopen(bin,'wb');
fwrite(fd,s,'float32')
fclose(fd);
