### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 1f1d6aef-78be-4bb5-901f-d352392b8f00
#Question #3
using LinearAlgebra, Random

# ╔═╡ 2c57d2ae-40c3-4685-8032-b7384eb25d9e
#Writing my Classical and Modified Gram-Schmid functions
function cgs(A)
	m, n = size(A)
	Q̂ = zeros(eltype(A), m, n)
	R̂ = zeros(eltype(A), n, n)
	for j = 1:n
		vj = A[:, j]
		for i = 1 : j - 1
			R̂[i, j] = Q̂[:, i]' * A[:, j]
			vj -= R̂[i, j] * Q̂[:, i] 
		end
		R̂[j, j] = norm(vj)
		Q̂[:, j] = vj / R̂[j, j]
	end
	return Q̂, R̂
end

# ╔═╡ cbdb1e4d-5e85-4c02-a4aa-6ae72128f463
function mgs(A)
    m, n = size(A)
    Q̂ = zeros(eltype(A), m, n)
    R̂ = zeros(eltype(A), n, n)
    V = copy(A)
    for j = 1:n
        R̂[j, j] = norm(V[:, j])
        Q̂[:, j] = V[:, j] / R̂[j, j]
        for i = j+1:n
            R̂[j, i] = dot(Q̂[:, j], V[:, i])
            V[:, i] -= R̂[j, i] * Q̂[:, j]
        end
    end
    return Q̂, R̂
end

# ╔═╡ 6384b06d-e698-49d6-9f6c-0596e9408f01
#here we are now going to test A1, A2, A3 for CGS and MGS
m = 100

# ╔═╡ 5014dab2-d206-47fe-96a8-096585bfe943
U = qr(randn(m)).Q

# ╔═╡ e09bf5d9-a2fa-4838-945d-8453efd7c79b
V = qr(randn(m)).Q

# ╔═╡ ee16f45a-d82a-4cef-8baa-f7b554b14f72
norm(V' * V - I)

# ╔═╡ c6a260d4-2316-49cf-a804-8f30a86318eb
norm(U' * U - I)

# ╔═╡ dd521fb9-5467-4808-aad5-61ad04f946d7
s1 = 10 .^ range(1, stop=-1, length=m)

# ╔═╡ 0037214f-424a-4f94-98eb-72df86d28769
s2 = 10 .^ range(1, stop=-5, length=m)

# ╔═╡ a3b43a4f-6021-41ba-aa8d-901f135d2308
s3 = 10 .^ range(1, stop=-25, length=m)

# ╔═╡ dde3f9f6-0d76-481e-9dd3-52d23142e1ed
A1_CGS = U * Diagonal(s1) * V'

# ╔═╡ 0d587928-6fe5-4c51-9eb6-40dd8b842caf
A2_CGS = U * Diagonal(s2) * V'

# ╔═╡ 690d8c9d-7ddf-4b1c-aeb1-492fb74a6243
A3_CGS = U * Diagonal(s3) * V'

# ╔═╡ 5966e6b4-2366-43b9-8e33-3853f8ba8983
#Now we will check the accuracy

# ╔═╡ 261fa5b1-26df-43f7-a1a6-cd4efc7f188b
Q̂_1cgs, R̂_1cgs = cgs(A1_CGS)

# ╔═╡ 8c1b8868-b291-4698-b6ec-bafa7ec02eeb
Q̂_2cgs, R̂_2cgs = cgs(A2_CGS)

# ╔═╡ a74463b1-cc94-4e76-b76f-db40df165e36
Q̂_3cgs, R̂_3cgs = cgs(A3_CGS)

# ╔═╡ cef05d9f-429a-4622-ad4d-2cd4b9605a02
norm(A1_CGS - Q̂_1cgs * R̂_1cgs )

# ╔═╡ 8fc36efe-7748-4eea-9af9-f5ed61f05428
norm(A2_CGS - Q̂_2cgs * R̂_2cgs )

# ╔═╡ df2fffda-0e8e-48bc-abc6-239670b54d25
norm(A3_CGS - Q̂_3cgs * R̂_3cgs )

# ╔═╡ dec28c99-de64-4968-892a-deb31568ec88
norm(Q̂_1cgs' * Q̂_1cgs - I)

# ╔═╡ ec0e1d21-0487-4b17-ad78-835315195bcb
norm(Q̂_2cgs' * Q̂_2cgs - I)

# ╔═╡ aff4f2dc-6d65-436c-9a16-71a24eb2e69a
norm(Q̂_3cgs' * Q̂_3cgs - I)

# ╔═╡ 0cd45f64-1104-4d2d-87a3-f86185640ab1
############################################################################

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
"""

# ╔═╡ Cell order:
# ╠═1f1d6aef-78be-4bb5-901f-d352392b8f00
# ╠═2c57d2ae-40c3-4685-8032-b7384eb25d9e
# ╠═cbdb1e4d-5e85-4c02-a4aa-6ae72128f463
# ╠═6384b06d-e698-49d6-9f6c-0596e9408f01
# ╠═5014dab2-d206-47fe-96a8-096585bfe943
# ╠═e09bf5d9-a2fa-4838-945d-8453efd7c79b
# ╠═ee16f45a-d82a-4cef-8baa-f7b554b14f72
# ╠═c6a260d4-2316-49cf-a804-8f30a86318eb
# ╠═dd521fb9-5467-4808-aad5-61ad04f946d7
# ╠═0037214f-424a-4f94-98eb-72df86d28769
# ╠═a3b43a4f-6021-41ba-aa8d-901f135d2308
# ╠═dde3f9f6-0d76-481e-9dd3-52d23142e1ed
# ╠═0d587928-6fe5-4c51-9eb6-40dd8b842caf
# ╠═690d8c9d-7ddf-4b1c-aeb1-492fb74a6243
# ╠═5966e6b4-2366-43b9-8e33-3853f8ba8983
# ╠═261fa5b1-26df-43f7-a1a6-cd4efc7f188b
# ╠═8c1b8868-b291-4698-b6ec-bafa7ec02eeb
# ╟─a74463b1-cc94-4e76-b76f-db40df165e36
# ╠═cef05d9f-429a-4622-ad4d-2cd4b9605a02
# ╠═8fc36efe-7748-4eea-9af9-f5ed61f05428
# ╠═df2fffda-0e8e-48bc-abc6-239670b54d25
# ╠═dec28c99-de64-4968-892a-deb31568ec88
# ╠═ec0e1d21-0487-4b17-ad78-835315195bcb
# ╠═aff4f2dc-6d65-436c-9a16-71a24eb2e69a
# ╠═0cd45f64-1104-4d2d-87a3-f86185640ab1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
