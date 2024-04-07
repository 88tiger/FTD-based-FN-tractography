function seeds = GenSeeds(lgnposition,FiberCount)

    multi = round(FiberCount/size(lgnposition,1));
    seeds_round = repmat(lgnposition,multi,1);
    seeds = seeds_round+rand(size(seeds_round))-[0.5,0.5,0.5];
end