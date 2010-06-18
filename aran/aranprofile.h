/* LIBARAN - Fast Multipole Method library
 * Copyright (C) 2006-2007 Pierre Gay
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __ARAN_PROFILE_H__
#define __ARAN_PROFILE_H__

#include <vsg/vsgd.h>
#include <aran/aran.h>
#include <aran/aransolver2d.h>
#include <aran/aransolver3d.h>

#include <aranrusage.h>
#include <aranpoly1d.h>

typedef void (*AranParticleInitFunc2d) (VsgPoint2 particle);
typedef gpointer (*AranDevelopmentNewFunc) (guint8 posdeg, guint8 negdeg);

gdouble aran_profile_p2p_2d (AranParticle2ParticleFunc2d p2p,
                             AranParticleInitFunc2d init,
                             VsgPoint2 p1,
                             VsgPoint2 p2);

gdouble aran_profile_nearfunc_2d (AranParticle2ParticleFunc2d p2p,
                                  AranParticleInitFunc2d init,
                                  VsgPoint2 p1,
                                  VsgPoint2 p2,
                                  gint maxbox);

gdouble aran_profile_p2m_2d (AranParticle2MultipoleFunc2d p2m,
                             AranZeroFunc init,
                             VsgPoint2 p,
                             VsgPRTree2dNodeInfo *nodeinfo,
                             gpointer dst);

gdouble aran_profile_m2m_2d (AranMultipole2MultipoleFunc2d m2m,
                             AranZeroFunc init,
                             VsgPRTree2dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree2dNodeInfo *dstinfo,
                             gpointer dst);


gdouble aran_profile_m2l_2d (AranMultipole2LocalFunc2d m2l,
                             AranZeroFunc init,
                             VsgPRTree2dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree2dNodeInfo *dstinfo,
                             gpointer dst);

gdouble aran_profile_l2l_2d (AranLocal2LocalFunc2d l2l,
                             AranZeroFunc init,
                             VsgPRTree2dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree2dNodeInfo *dstinfo,
                             gpointer dst);

gdouble aran_profile_l2p_2d (AranLocal2ParticleFunc2d l2p,
                             AranParticleInitFunc2d init,
                             VsgPRTree2dNodeInfo *nodeinfo,
                             gpointer dst,
                             VsgPoint2 p);

gdouble aran_poly1d_profile_p2m_2d (AranParticle2MultipoleFunc2d p2m,
                                    AranZeroFunc init,
                                    VsgPoint2 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_p2m_2d_samples (AranParticle2MultipoleFunc2d p2m,
                                            AranZeroFunc init,
                                            VsgPoint2 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

gdouble aran_poly1d_profile_m2m_2d (AranMultipole2MultipoleFunc2d m2m,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_m2m_2d_samples (AranMultipole2MultipoleFunc2d m2m,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

gdouble aran_poly1d_profile_m2l_2d (AranMultipole2LocalFunc2d m2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_m2l_2d_samples (AranMultipole2LocalFunc2d m2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);


gdouble aran_poly1d_profile_l2l_2d (AranLocal2LocalFunc2d l2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_l2l_2d_samples (AranLocal2LocalFunc2d l2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

gdouble aran_poly1d_profile_l2p_2d (AranLocal2ParticleFunc2d l2p,
                                    AranParticleInitFunc2d init,
                                    VsgPoint2 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_l2p_2d_samples (AranLocal2ParticleFunc2d l2p,
                                            AranParticleInitFunc2d init,
                                            VsgPoint2 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

typedef void (*AranParticleInitFunc3d) (VsgPoint3 particle);

gdouble aran_profile_p2p_3d (AranParticle2ParticleFunc3d p2p,
                             AranParticleInitFunc3d init,
                             VsgPoint3 p1,
                             VsgPoint3 p2);

gdouble aran_profile_nearfunc_3d (AranParticle2ParticleFunc3d p2p,
                                  AranParticleInitFunc3d init,
                                  VsgPoint3 p1,
                                  VsgPoint3 p2,
                                  gint maxbox);

gdouble aran_profile_p2m_3d (AranParticle2MultipoleFunc3d p2m,
                             AranZeroFunc init,
                             VsgPoint3 p,
                             VsgPRTree3dNodeInfo *nodeinfo,
                             gpointer dst);

gdouble aran_profile_m2m_3d (AranMultipole2MultipoleFunc3d m2m,
                             AranZeroFunc init,
                             VsgPRTree3dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree3dNodeInfo *dstinfo,
                             gpointer dst);


gdouble aran_profile_m2l_3d (AranMultipole2LocalFunc3d m2l,
                             AranZeroFunc init,
                             VsgPRTree3dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree3dNodeInfo *dstinfo,
                             gpointer dst);

gdouble aran_profile_l2l_3d (AranLocal2LocalFunc3d l2l,
                             AranZeroFunc init,
                             VsgPRTree3dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree3dNodeInfo *dstinfo,
                             gpointer dst);

gdouble aran_profile_l2p_3d (AranLocal2ParticleFunc3d l2p,
                             AranParticleInitFunc3d init,
                             VsgPRTree3dNodeInfo *nodeinfo,
                             gpointer dst,
                             VsgPoint3 p);

gdouble aran_poly1d_profile_p2m_3d (AranParticle2MultipoleFunc3d p2m,
                                    AranZeroFunc init,
                                    VsgPoint3 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_p2m_3d_samples (AranParticle2MultipoleFunc3d p2m,
                                            AranZeroFunc init,
                                            VsgPoint3 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

gdouble aran_poly1d_profile_m2m_3d (AranMultipole2MultipoleFunc3d m2m,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_m2m_3d_samples (AranMultipole2MultipoleFunc3d m2m,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

gdouble aran_poly1d_profile_m2l_3d (AranMultipole2LocalFunc3d m2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_m2l_3d_samples (AranMultipole2LocalFunc3d m2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

gdouble aran_poly1d_profile_l2l_3d (AranLocal2LocalFunc3d l2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_l2l_3d_samples (AranLocal2LocalFunc3d l2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

gdouble aran_poly1d_profile_l2p_3d (AranLocal2ParticleFunc3d l2p,
                                    AranParticleInitFunc3d init,
                                    VsgPoint3 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples);

gdouble aran_poly1d_profile_l2p_3d_samples (AranLocal2ParticleFunc3d l2p,
                                            AranParticleInitFunc3d init,
                                            VsgPoint3 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples);

#endif /* __ARAN_PROFILE_H__ */
