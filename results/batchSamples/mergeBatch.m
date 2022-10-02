% Merge the batches of samples for the reference computation

merged_o = [];
merged_Y = [];
level = 4;

for i=0:99
    file_name = sprintf('batchSamples/batchSamples-%d.mat', i);
    load(file_name);
    merged_o = cat(1, merged_o, o);
    merged_Y = cat(1, merged_Y, Y);
end

save 