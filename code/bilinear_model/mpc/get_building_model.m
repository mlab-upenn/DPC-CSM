function [bm] = get_building_model(build)

    
cd(['InData/',build.sys,'_occup'])
    multiZapprox = 0;
    filename = sprintf('%s-Or_%s_%s_%s_%s_%s_d_multiZapprox%s.txt', build.sys, build.dir, build.type, build.env, build.win, build.int, num2str(multiZapprox));
    building = filename;
    if strcmp(build.sys, 'e01')   % always have multiZapprox == 0 (implicitly)
        range = 'A11..L22';
        bm.A = dlmread(building, '\t', range);
        range = 'A25..F36';    
        bm.Bu = dlmread(building, '\t', range);
        range = 'A39..H50';
        bm.Bv = dlmread(building, '\t', range);
        range = 'A53..H64';
        bm.Bvu(:,:,1) = dlmread(building, '\t', range);
        range = 'A67..H78';
        bm.Bvu(:,:,2) = dlmread(building, '\t', range);
        range = 'A81..H92';
        bm.Bvu(:,:,3) = dlmread(building, '\t', range);
        range = 'A95..H106';
        bm.Bvu(:,:,4) = dlmread(building, '\t', range);
        range = 'A109..H120';
        bm.Bvu(:,:,5) = dlmread(building, '\t', range);
        range = 'A123..H134';
        bm.Bvu(:,:,6) = dlmread(building, '\t', range);
        range = 'A137..L148';
        bm.Bxu(:,:,1) = dlmread(building, '\t', range);
        range = 'A151..L162';
        bm.Bxu(:,:,2) = dlmread(building, '\t', range);
        range = 'A165..L176';
        bm.Bxu(:,:,3) = dlmread(building, '\t', range);
        range = 'A179..L190';
        bm.Bxu(:,:,4) = dlmread(building, '\t', range);
        range = 'A193..L204';
        bm.Bxu(:,:,5) = dlmread(building, '\t', range);
        range = 'A207..L218';
        bm.Bxu(:,:,6) = dlmread(building, '\t', range);
        range = 'A221..L223';
        bm.C = dlmread(building, '\t', range);
        range = 'A226..F228';
        bm.Du = dlmread(building, '\t', range);
        range = 'A231..H233';
        bm.Dv = dlmread(building, '\t', range);
        range = 'A236..H238';
        bm.Dvu(:,:,1) = dlmread(building, '\t', range);
        range = 'A241..H243';
        bm.Dvu(:,:,2) = dlmread(building, '\t', range);
        range = 'A246..H248';
        bm.Dvu(:,:,3) = dlmread(building, '\t', range);
        range = 'A251..H253';
        bm.Dvu(:,:,4) = dlmread(building, '\t', range);
        range = 'A256..H258';
        bm.Dvu(:,:,5) = dlmread(building, '\t', range);
        range = 'A261..H263';
        bm.Dvu(:,:,6) = dlmread(building, '\t', range);
    end
    cd ..
    cd ..
    
    % remove u4
    bm.Bu(:,4) = [];
    bm.Bxu(:,:,4) = [];
    bm.Bvu(:,:,4) = [];
    bm.Du(:,4) = [];
    bm.Dvu(:,:,4) = [];
    
    % remove u3
    bm.Bu(:,3) = [];
    bm.Bxu(:,:,3) = [];
    bm.Bvu(:,:,3) = [];
    bm.Du(:,3) = [];
    bm.Dvu(:,:,3) = [];
end