# -*- coding: utf-8 -*-
'''
Author @55-AA
19 July, 2016
'''
import copy

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def mod_inv(a, m):
    g, x, y = egcd(a, m)
    assert g == 1
    return x % m

class GaussMatrix:
    def __init__(self, matrix, mod):
        self.matrix = copy.deepcopy(matrix)
        self.d = None

        self.r = len(matrix)
        self.c = len(matrix[0])
        self.N = len(matrix[0]) - 1
        self.mod = mod
        self.count = 1
        self.error_str = "unknown error"
        
    def swap_row(self, ra, rb):
        (self.d[ra], self.d[rb]) = (self.d[rb], self.d[ra])

    def swap_col(self, ca, cb):
        for j in range(self.r):
            (self.d[j][ca], self.d[j][cb]) = (self.d[j][cb], self.d[j][ca])

    def inv_result(self, r, n):
        b = self.d[n][self.N]
        a = self.d[n][n]
        m = self.mod
        k = gcd(a, m)            
        for j in xrange(n + 1, self.N):
            b = (b - (self.d[n][j] * r[j] % m)) % m

        if 1 == k:
            return [mod_inv(a, m) * b % m]
        else:
            if k == gcd(k, b):
                a /= k
                b /= k
                m /= k

                x0 = mod_inv(a, m) * b % m
                x = []
                for i in xrange(k):
                    x.append(x0 + m*i)
                return x
        return None

    def find_min_gcd_row_col(self, i, j):
        for k in xrange(i, self.r):
            for l in xrange(j, self.c - 1):
                if(1 == gcd(self.d[k][l], self.mod)):
                    return [k, l]


        def add_min_gcd(a, b, m):
            r = [m, 1]
            g = gcd(a, b)
            if g:
                i = a / g
                for j in xrange(i):
                    g = gcd((a + j * b) % m, m)
                    if g < r[0]:
                        r[0] = g
                        r[1] = j
                    if g == 1:
                        break
            return r

        r = [self.mod, 1, i, i + 1, j]
        for k in xrange(i, self.r):
            for kk in xrange(k+1, self.r):
                for l in range(j, self.c - 1):
                    rr = add_min_gcd(self.d[k][l], self.d[kk][l], self.mod)
                    if rr[0] < r[0]:
                        r[0] = rr[0]
                        r[1] = rr[1]
                        r[2] = k
                        r[3] = kk
                        r[4] = l
                        pass
                    if(1 == rr[0]):
                        break
        g = r[0]
        n = r[1]
        k = r[2]
        kk = r[3]
        l = r[4]

        if n and g < self.mod:
            self.d[k] = map(lambda x, y : (x + n*y)%self.mod, self.d[k], self.d[kk])
        return [k, l]
        
    def mul_row(self, i, k, j):
        a = self.d[k][j]
        b = self.d[i][j]

        def get_mul(a, b, m):
            k = gcd(a, m)
            if 1 == k:
                return mod_inv(a, m) * b % m
            else:
                if k == gcd(k, b):
                    return mod_inv(a/k, m/k) * (b/k) % (m/k)
            return None

        if b:
            mul = get_mul(a, b, self.mod)
            if None == mul:
                print_matrix(self.d)
                assert(mul != None)
            self.d[i] = map(lambda x, y : (y - x*mul) % self.mod, self.d[k], self.d[i])


    def gauss(self):
        self.d = copy.deepcopy(self.matrix)
        for i in xrange(self.r):
            for j in xrange(self.c):
                self.d[i][j] = self.matrix[i][j] % self.mod

        if self.r < self.N:
            self.d.extend([[0]*self.c]*(self.N - self.r))          

        index = [x for x in xrange(self.N)]
        for i in range(self.N):
            tmp = self.find_min_gcd_row_col(i, i)
            if(tmp):
                self.swap_row(i, tmp[0])
                (index[i], index[tmp[1]]) = (index[tmp[1]], index[i])
                self.swap_col(i, tmp[1])
            else:
                self.error_str = "no min"
                return None

            for k in range(i + 1, self.r):
                self.mul_row(k, i, i)

        if self.r > self.N:
            for i in xrange(self.N, self.r):
                for j in xrange(self.c):
                    if self.d[i][j]:
                        self.error_str = "r(A) != r(A~)"
                        return None

        for i in xrange(self.N):
            self.count *= gcd(self.d[i][i], self.mod)

        if self.count > 100:
            self.error_str = "solution too more:%d" % (self.count)
            return None            

        result = [[None]*self.N]
        for i in range(self.N - 1, -1, -1):
            new_result = []
            for r in result:
                ret = self.inv_result(r, i)
                if ret:
                    for rr in ret:
                        l = r[:]
                        l[i] = rr
                        new_result.append(l)

                else:
                    self.error_str = "no inv:i=%d" % (i)
                    return None

            result = new_result

        for i in xrange(len(result)) :
            def xchg(a, b):
                result[i][b] = a
            map(xchg, result[i][:], index)

        return result

def run_test(mod, matrix):
    g = GaussMatrix(matrix, mod)
    ret = g.gauss()
    if not ret:
        print "error:"
        print_matrix(g.d)
        print "error_str:", g.error_str
    else:
        print ret[0]

if __name__ == "__main__":
    mod = 6
    matrix = [
    [1, 0, 2, 6 - 0],
    [2, 0, 0, 6 - 2],
    [0, 4, 3, 6 - 3]
    ]
    run_test(mod, matrix)
