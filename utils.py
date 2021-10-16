import random

class BaseQualities:
    @staticmethod
    def QToPhred33(Q: int):
        return chr(Q + 33)

    @staticmethod
    def phred33ToQ(qual: int):
        return ord(qual) + 33


def random_embedded_genome(p: str, l: int = 100):
    rand_pos = random.randint(0, l - len(p))
    seq = ''
    for i in range(l):
        if i == rand_pos:
            seq += p
        else:
            seq += random.choice('ATCG')
    return seq, rand_pos


class StringMatching:
    @staticmethod
    def naive(p: str, t: str):
        l = 0
        r = len(p)
        for i in range(len(t) - len(p) + 1):
            match = True
            for j in range(len(p)):
                if t[l + j] == p[j]:
                    # Same chars
                    continue
                else:
                    match = False
                    break
            if match:
                return {
                    'match_found': True,
                    'start_location': l
                }
            else:
                l += 1
                r += 1
        return {
                'match_found': False,
                'start_location': -1
            }

