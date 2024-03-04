function[reg_stack, tform] = corrRegMax(stackA, stackB)


max_z_A = max(max(stackA,[],4),[],3);
max_z_B = max(max(stackB,[],4),[],3);

tform = imregcorr(max_z_B,max_z_A)

Rfixed = imref2d(size(max_z_A));

sz = size(stackB);

reg_stack = zeros(size(stackA));

for channels = 1:sz(4)
    for z = 1:sz(3)
        
    movingimg = stackB(:,:,z,channels);

    movingReg = imwarp(movingimg,tform,'OutputView',Rfixed);
    
    reg_stack(:,:, z, channels) = movingReg; 

    end
end

end