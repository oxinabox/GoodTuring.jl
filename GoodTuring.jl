#' This is an implementation of Simple Good Turing Smothing,
#' as described by 


module GoodTuring

	using DataFrames, StatsBase#, GLM
	export simpleGoodTuring


	function simpleGoodTuring(speciesCountDict::Dict)
	    speciesCountVec = collect(values(speciesCountDict))
	        
	    totalCounts = sum(speciesCountVec)
	    cofcDict = countmap(speciesCountVec)
	    r = sort(collect(keys(cofcDict)))
	
	    N = size(r,1)
	    Nr = [cofcDict[r[i]] for i in 1:N]
	
	    p0 = haskey(cofcDict, 1.0) ? cofcDict[1.0] / totalCounts : 0.0
	    
	    Z = sgtZ(r,Nr)
	    logr = map(log,r)
	    logZ = map(log,Z)
	
	    X = hcat(ones(N), logr)
	    Y = copy(logZ )
	    coefs = X\Y
	    intercept = coefs[1]
	    slope = coefs[2]
	
	    useY = false
	    rSmooth = Array{Float64}(N)
	    for i in 1:N
	        @inbounds thisr = r[i]
	        
	        #y = ((thisr+1.0)^(slope+1.0))/(thisr^slope)
	        #The above is the much simplified form of the below (Performance identical output differs by 10^-16)
	        y = (thisr+1.0) * exp(slope * log(thisr+1.0) + intercept) / exp(slope * log(thisr) + intercept)
	
	        if !in(thisr+1, r)
	            useY = true
	        end
	
	        if useY
	            rSmooth[i] = y
	        else
	            x = (thisr+1) * cofcDict[thisr + 1]/cofcDict[thisr]
	            thisNr = cofcDict[thisr]
	            thisNr1 = cofcDict[thisr+1]
	
	            t = 1.96 * ((thisr+1)^2) * (thisNr1 / thisNr^2) * (1 + (thisNr1 / thisNr))
	
	            if abs(x-y) > t
	                @inbounds rSmooth[i] = x
	            else
	                useY = true
	                @inbounds rSmooth[i] = y
	            end
	        end
	    end
	
	    smoothTot = sum(Nr.*rSmooth)
	    sgtProb  = (1.0 - p0) .* (rSmooth/smoothTot)
	    sgtProbDict = Dict([r[i] => sgtProb[i] for i in 1:N])
	    sgtDict = Dict([sp=>sgtProbDict[speciesCountDict[sp]] for sp in keys(speciesCountDict)])
	
	    sgtDict, sgtProbDict, p0
	end

	function simpleGoodTuring(x::DataFrames.DataArray)
		speciesCountDict = countmap(x)
		df, p0 = simpleGoodTuring(speciesCountDict)
		return df, p0
	end

	function simpleGoodTuring(x::Array)
		speciesCountDict = countmap(x)
		df, p0 = simpleGoodTuring(speciesCountDict)
		return df, p0
	end


	function sgtZ(r::Array, Nr::Array)
	    j = r
	    i = [0, j[1:end-1]]
	    lastK = 2*j[end] - i[end]
	    k = [j[2:end], lastK]
	    Float64[(2*Nr[iter])/(k[iter]-i[iter]) for iter = 1:length(j)]
	end



end
