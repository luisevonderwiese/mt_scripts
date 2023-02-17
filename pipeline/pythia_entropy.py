    def entropy_matrix(self) -> float:

        def _remove_gaps_from_sequence(seq):
            for char in GAP_CHARS:
                seq = seq.replace(char, "")
            return seq

        entropies = []
        flat_matrix = []
        for i in range(len(self.msa)):
            flat_matrix += _remove_gaps_from_sequence(self.msa[i, :]).upper()
        entropy = 0

        for char in STATE_CHARS:
            count = str.count(flat_matrix, char)
            if count == 0:
                entropy_x = 0
            else:
                prob = count / len(flat_matrix)
                entropy_x = prob * math.log2(prob)

            entropy += entropy_x

        entropy = -entropy

        assert (
                entropy >= 0
        ), f"Entropy negative, check computation. Entropy is {entropy}"

        return entropy


    def row_entropies(self) -> List[float]:

        def _remove_gaps_from_sequence(seq):
            for char in GAP_CHARS:
                seq = seq.replace(char, "")
            return seq

        entropies = []
        for i in range(len(self.msa)):
            row = _remove_gaps_from_sequence(self.msa[i, :]).upper()
            entropy = 0

            for char in STATE_CHARS:
                count = str.count(row, char)
                if count == 0:
                    entropy_x = 0
                else:
                    prob = count / len(row)
                    entropy_x = prob * math.log2(prob)

                entropy += entropy_x

            entropy = -entropy

            assert (
                    entropy >= 0
            ), f"Entropy negative, check computation. Entropy is {entropy}"

            entropies.append(entropy)
        return entropies

    def column_entropies(self) -> List[float]:
        """Returns the shannon entropy (in bits) for each site in the MSA.

        Returns:
            column_entropies (List[float]): List of shannon entropies corresponding to the MSA site. Each per-site entropy is >= 0.
        """

        def _remove_gaps_from_sequence(seq):
            for char in GAP_CHARS:
                seq = seq.replace(char, "")
            return seq

        entropies = []
        for i in range(self.msa.get_alignment_length()):
            column = _remove_gaps_from_sequence(self.msa[:, i]).upper()
            entropy = 0

            for char in STATE_CHARS:
                count = str.count(column, char)
                if count == 0:
                    entropy_x = 0
                else:
                    prob = count / len(column)
                    entropy_x = prob * math.log2(prob)

                entropy += entropy_x

            entropy = -entropy

            assert (
                    entropy >= 0
            ), f"Entropy negative, check computation. Entropy is {entropy}"

            entropies.append(entropy)
        return entropies

    def entropy(self) -> float:
        """Returns the shannon entropy (in bits) of the MSA.

        Returns:
            entropy (float): Shannon entropy of the MSA. The entropy is >= 0.
        """
        return statistics.mean(self.column_entropies())



    def entropy_row(self) -> float:
        """Returns the shannon entropy (in bits) of the MSA.

        Returns:
            entropy (float): Shannon entropy of the MSA. The entropy is >= 0.
        """
        return statistics.mean(self.row_entropies())

    def bollback_multinomial(self) -> float:
        """Returns the bollback multinomial statistic of the MSA.

        According to Bollback, JP: Bayesian model adequacy and choice in phylogenetics (2002)

        Returns:
            bollback (float): The bollback multionomial statistic of the MSA. The bollback multinomial statistic is <= 0.
        """
        msa_length = self.number_of_sites()

        sites = []
        for i in range(msa_length):
            sites.append(self.msa[:, i])

        site_counts = Counter(sites)
        mult = 0
        for i in site_counts:
            N_i = site_counts[i]
            mult += N_i * math.log(N_i)

        mult = mult - msa_length * math.log(msa_length)
        return mult


    def pattern_entropy(self) -> float:
        msa_length = self.number_of_sites()

        sites = []
        for i in range(msa_length):
            sites.append(self.msa[:, i])

        site_counts = Counter(sites)
        entropy = 0
        for site in site_counts:
            count = site_counts[site]
            if count == 0:
                entropy_x = 0
            else:
                prob = count / len(site_counts)
                entropy_x = prob * math.log2(prob)

            entropy += entropy_x

        entropy = -entropy
        assert (
                entropy >= 0
        ), f"Entropy negative, check computation. Entropy is {entropy}"

        return entropy


