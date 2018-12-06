"""
`coinflip'
===========
Generate a random graph by coin flipping.  The input is a 
-----
- `p` : square matrix with entries between 0 and 1. 
Functions
---------
- `coinflip(p::Float64)` 
Output
------
- A 0 or 1 adjacency matrix for the random graph.


"""

function coinflip(p::Float64)
	n = size(p,1)
	type = eltype(p)
	n == size(p,2) ? nothing : throw(AssertionError("Matrix not square"))   # make sure we have a square input
	A = zeros(type, n, n)                                                   # create an empty adjacency matrix
	for i=1:n
		for j=1:n                                                           # fill in each entry as 1 with prob p[i,j] 
			A[i,j] = rand() < p[i,j]
		end
	end
	return A
end

#println(coinflip(.25*ones(8,8)))  # generate a random Erdos-Renyi graph

"""
`ball_drop_er`
===============
Generate a random Erdos-Renyi graph by ball-dropping.
Input
-----
- `n` : the number of nodes
- `p` : the probablity of an edge
Functions
---------
- `ball_drop_er`(n::Int64,p::Float64)
Output
-------
- List of directed edges.

"""

function ball_drop_er(n::Int64,p::Float64)
	m = rand(Binomial(n*n,p))              # the number of edges
	edges = Set{typeof((n,n))}()           # create a set of type n
	while length(edges) < m
		push!(edges,(rand(0:n-1), rand(0:n-1)))  # use indices in 0, n-1
	end
	return edges
end

#println(ball_drop_er(8,.25))  # 8 node Erdos-Renyi with probablity 0.25

"""
`grass_hop_er`
===============
Generate a random Erdos-Renyi graph by grass-hopping. 
Input
-----
- `n` : the number of nodes
- `p` : the probablity of an edge
Functions
---------
- `grass_hop_er(n::Int64,p::Float64)`
Output
------
- The result is a list of directed edges.

"""

function grass_hop_er(n::Int64,p::Float64)
	edgeindex = -1                 # we label edges from 0 to n^2-1
	gap = rand(Geometric(p))       # first distance to edge
	edges = Set{typeof((n,n))}()   # make a set of type n
	while edgeindex+gap < n*n      # check to make sure we have a valid index
		edgeindex += gap             # increment the index
		src = floor(edgeindex/n)     # use floor divison, gives src in [0, n-1]
		dst = edgeindex - (n*src)    # identify the column
		push!(edges,(src,dst))
		gap = rand(Geometric(p))     # generate the next gap
	end
	return edges
end

#println(grass_hop_er(8,0.25))

"""
`sbm2`
======
Generate edges for a stochastic block model with two blocks.

Input
-----
- `n1` : Size of block
- `n2` : Size of block
- `p`  : within group probablity
- `q`  : between group probablity
Functions
---------
- `sbm2(n1::Int64,n2::Int64,p::Float64,q::Float64)`

Output
------
- Edges for a stochastic block model with two blocks

"""
function sbm2(n1::Int64,n2::Int64,p::Float64,q::Float64)
	edges = grass_hop_er(n1,p)               #genertate the n1-by-n1 block
	for (i,j) in grass_hop_er(n2,p)          #n2-by-n2 block
		push!(edges,(i+n1,j+n1))
	end
	for (i,j) in grass_hop_er(max(n1,n2),q)
		if i<n1 && j<n2
	push!(edges,(i,j+n1))
		end
	end
	for (i,j) in grass_hop_er(max(n1,n2),q) 
		if i<n2 && j<n1
			 push!(edges,(i+n1,j))
		end
	end
	return edges
end  

#println(sbm2(20,15,.5,.1)) 


"""
- `unrank`
===========
Input
-----
- `C`: a multiset represented by a list
- `n`: the lexicographic rank to find
Functions
---------
- `ndseq_to_counter(seq)`
- `counter_to_ndseq(mset, keys)`
- `num_multiset_permutations(mset)`
- `unrank_mset_counter(mset,keys,n)`
- `unrank(seq,n)`

Output
------
The nth permutation of C in lexicographic order
"""

function ndseq_to_counter(seq)
	mset = Dict{Int64, Int64}()
	for c in seq
		#get a value with a default
		#of zero if it isn't there
		mset[c] = get(mset, c, 0) +1
	end
	return mset, sort(collect(keys(mset)))
end

function counter_to_ndseq(mset, keys)
	seq = Int64[]
		for k in keys #keys in sorted order
					  #append k mset[k] times
			for v in range(1, stop = mset[k])
				append!(seq, k)
			end
		end
	return seq
end

function num_multiset_permutations(mset)
	count = factorial(sum(values(mset)))
	return count
end

function unrank_mset_counter(mset,keys,n)
	if n == 0		#easy if rank == 0
		return counter_to_ndseq(mset,keys)
	end
	for s in keys    	#else find prefix key
		mset[s] -= 1	#decrease count of s
						# determine number of prefixes with s
		place = num_multiset_permutations(mset)
		if place > n
								#then the first element is s
			if mset[s] == 0		#remove the key if count is 0
				keys = setdiff(keys, s)
			end
			suffix = unrank_mset_counter(mset, keys, n) #recurse
			pushfirst!(suffix, s)			#append s
			return suffix
		else				#here it does not start with s
			mset[s] += 1	#restore the count
			n -= place		#update search offset
		end
	end
	throw(ValueError("rank too large"))		
end

function unrank(seq,n)
	mset, keys = ndseq_to_counter(seq)
	return unrank_mset_counter(mset,keys,n)
end

println(unrank([0,1,1,3], 0))
println(unrank([0,1,1,3], 1))
println(unrank([0,1,1,3], 2))


"""
`map_mult_to_kron`
===================
Map a multi-index from the mult. table 
table to a row and column in the Kronecker 
matrix.
Input
-----	
- `mind` : the multi-index for the mult table
- `n`    : the size of the initiator matrix K
Functions
---------
- `map_mult_to_kron(mind,n::Int64)`
- `multiindex_to_linear(mind, n2::Int64)`
- `morton_decode(I::Int64,n::Int64)`
Output
------


"""

function map_mult_to_kron(mind,n::Int64)
	I = multiindex_to_linear(mind, n*n)
	return morton_decode(I,n)
end

function multiindex_to_linear(mind, n2::Int64)
	I = 0
	base = 1
	for i in range(length(mind), step = -1, stop = 1)
		I +=mind[i]*base
		base *= n2
	end
	return I
end

function morton_decode(I::Int64,n::Int64)
	row = 0
	rowbase = 1
	col = 0
	colbase = 1
	i = 0
	while I > 0
		digit = I%n
		I = fld(I, n)
		if i%2 == 0
			row += rowbase * digit
			rowbase *= n
		else
			col += colbase * digit
			colbase *= n
		end
		i += 1
	end
	return row, col
end

#println(map_mult_to_kron([1,3],2))
#println(map_mult_to_kron([4,0,7],3))


