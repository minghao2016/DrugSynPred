function [ stitchchD2D ] = Construct_stitchD2D()

    stitch = load('input/STITCH/drug_comb_mapped');
    value = stitch(:,3);
    value = value / norm(value);
    stitchchD2D = sparse(stitch(:,1),stitch(:,2),value,119,119);
    stitchchD2D = full(stitchchD2D);
end