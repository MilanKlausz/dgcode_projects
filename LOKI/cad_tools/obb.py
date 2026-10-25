"""Oriented boxes: separating-axis overlap test, containment, conversion between frames."""
import numpy as np


class Box:
    """centre c (3), axes A (3x3, rows = unit axes), half sizes h (3)"""
    def __init__(self, c, A, h, name='', mat='', meta=None):
        self.c = np.asarray(c, float); self.A = np.asarray(A, float); self.h = np.asarray(h, float)
        self.name = name; self.mat = mat; self.meta = meta or {}

    def corners(self):
        s = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)], float)
        return self.c + (s * self.h) @ self.A

    def transformed(self, R, t):
        """world = R @ local + t"""
        return Box(R @ self.c + t, (R @ self.A.T).T, self.h, self.name, self.mat, dict(self.meta))

    def contains(self, other, margin=0.0):
        P = other.corners() - self.c
        loc = P @ self.A.T
        return bool(np.all(np.abs(loc) <= self.h - margin + 1e-9))


def overlap_depth(a, b):
    """separating axis test of two oriented boxes: the penetration depth along the best separating direction
    (min over the 15 SAT axes of ra + rb - |d.L|); < 0: separated, by at least -depth"""
    axes = list(a.A) + list(b.A)
    for i in range(3):
        for j in range(3):
            c = np.cross(a.A[i], b.A[j])
            if np.linalg.norm(c) > 1e-9:
                axes.append(c / np.linalg.norm(c))
    d = b.c - a.c
    return min(np.sum(a.h * np.abs(a.A @ L)) + np.sum(b.h * np.abs(b.A @ L)) - abs(d @ L) for L in axes)


def separating_axis(a, b):
    """the SAT axis with the smallest penetration, oriented from a to b, and that penetration"""
    axes = list(a.A) + list(b.A)
    for i in range(3):
        for j in range(3):
            c = np.cross(a.A[i], b.A[j])
            if np.linalg.norm(c) > 1e-9:
                axes.append(c / np.linalg.norm(c))
    d = b.c - a.c
    pen = [np.sum(a.h * np.abs(a.A @ L)) + np.sum(b.h * np.abs(b.A @ L)) - abs(d @ L) for L in axes]
    k = int(np.argmin(pen)); L = axes[k]
    return (L if d @ L >= 0 else -L), pen[k]
