function[reg_stack, tform] = corrRegMax_morph(stackA, stackB)


max_z_A = max(stackA,[],3);
max_z_B = max(stackB,[],3);

tform = imregcorr(max_z_B,max_z_A)

Rfixed = imref2d(size(max_z_A));

sz = size(stackB);

reg_stack = zeros(size(stackA));

for z = 1:sz(3)
        
    movingimg = stackB(:,:,z);

    movingReg = imwarp(movingimg,tform,'OutputView',Rfixed);
    
    reg_stack(:,:, z) = movingReg; 

end

end