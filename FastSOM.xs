#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <math.h>
#include "FastSOM.h"
#include "proto.h"

U32 mysize() {
	U32 size = 0;
	char buf[30];
	snprintf(buf,30,"/proc/self/statm");
	FILE* f = fopen(buf,"r");
	if (f) fscanf(f,"%u",&size);
	fclose(f);
	return size;
}


NV _vector_distance(AV* V1, AV* V2) {
	NV	diff,sum;
	I32	w_ptr;

	sum = 0;

	for ( w_ptr=av_len(V2) ; w_ptr>=0 ; w_ptr-- ) {
		diff = SvNV(*av_fetch(V1, w_ptr, FALSE))
			- SvNV(*av_fetch(V2, w_ptr, FALSE));
		sum += diff * diff;
	}
	return sqrt(sum);
}

void _bmu_guts(SOM_GENERIC *som,AV *sample,IV *bx,IV *by,NV *bd) {
	IV		x,y,z,X,Y,Z;
	NV		sum,diff,distance;
	SOM_Map		*map;
	SOM_Array	*array;
	SOM_Vector	*vector;

	map = som->map;
	X = som->X;
	Y = som->Y;
	Z = som->Z;

	*bx = -1;
	*by = 0;
	*bd = 0.0;

	for ( x=0 ; x<X ; x++ ) {
		array = (SOM_Array*)(&map->array)[x];
		for ( y=0 ; y<Y ; y++ ) {
			vector = (SOM_Vector*)(&array->vector)[y];

			sum = 0;
			for ( z=0 ; z<Z ; z++ ) {
				diff = SvNV(*av_fetch(sample,z,0))
					- (&vector->element)[z];
				sum += diff * diff;
			}

			distance = sqrt(sum);

			if ( *bx < 0 )
				{ *bx = 0; *by = 0; *bd = distance; }
			if ( distance < *bd )
				{ *bx = x; *by = y; *bd = distance; }
		}
	}
}


/* http://www.ai-junkie.com/ann/som/som4.html */
void _adjust(SV* self,NV l,NV sigma,AV* unit,AV* v) {
	IV		x,y;
	I32		z,Z;
	NV		d,theta,vold,wold;
	MAGIC		*mg;
	SOM_Map		*map;
	SOM_Array	*array;
	SOM_Vector	*vector;
	SOM_GENERIC	*som;

	x = SvIV(*av_fetch(unit, 0, FALSE));
	y = SvIV(*av_fetch(unit, 1, FALSE));
	d = SvNV(*av_fetch(unit, 2, FALSE));
	theta = exp( -d*d/2/sigma/sigma );

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!\n");
	som = self2somptr(self,mg);

	map = som->map;
	array = (SOM_Array*)(&map->array)[x];
	vector = (SOM_Vector*)(&array->vector)[y];

	/* hmm.. casting IV to I32.. is that sane? */
	Z = (I32)som->Z;

	for ( z=0 ; z<Z ; z++ ) {
		wold = (&vector->element)[z];
		vold = SvNV(*av_fetch(v,z,FALSE));
		(&vector->element)[z] = (vold - wold) * l * theta + wold;
	}
}

void _adjustn(SOM_GENERIC* som,NV l,NV sigma,NV* n,AV* v) {
	IV		x,y,X,Y;
	I32		z,Z;
	NV		d,theta,vold,wold;
	SOM_Map		*map;
	SOM_Array	*array;
	SOM_Vector	*vector;

	map = som->map;
	X = som->X;
	Y = som->Y;

	for ( x=0 ; x<X ; x++ ) {
		array = (SOM_Array*)(&map->array)[x];
		for ( y=0 ; y<Y ; y++ ) {
			d = n[x*X+y];
			if (d < 0) continue;
			theta = exp( -d*d/2/sigma/sigma );
			vector = (SOM_Vector*)(&array->vector)[y];

			/* hmm.. casting IV to I32.. is that sane? */
			Z = (I32)som->Z;

			for ( z=0 ; z<Z ; z++ ) {
				wold = (&vector->element)[z];
				vold = SvNV(*av_fetch(v,z,FALSE));
				(&vector->element)[z] =
					(vold - wold) * l * theta + wold;
			}
		}
	}
}

void _som_new(const char* self,...) { croak("Dont use this class directly\n"); }
void _som_neighbors(SV* self,...)   { croak("unsupported"); }
void _som_as_string(SV* self,...)   { croak("unsupported"); }
void _som_as_data(SV* self,...)     { croak("unsupported"); }

void _som_bmu(SV* self, AV* sample) {
	IV		cx,cy;
	NV		cd;
	MAGIC		*mg;
	SOM_GENERIC	*som;
	dXSARGS;

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!\n");
	som = self2somptr(self,mg);

	_bmu_guts(som,sample,&cx,&cy,&cd);

	sp = mark;
	XPUSHs(sv_2mortal(newSViv(cx)));
	XPUSHs(sv_2mortal(newSViv(cy)));
	XPUSHs(sv_2mortal(newSVnv(cd)));
	PUTBACK;
}

SV* map(SV* self) {
        MAGIC		*mg;
	SOM_GENERIC	*som;

	if ( !(mg = selfmagic(self)) )
                croak("self has no magic!\n");
	som = self2somptr(self,mg);

	SvREFCNT_inc(som->map->ref);
	return som->map->ref;
}

SV* output_dim(SV* self) {
	MAGIC		*mg;
	SOM_GENERIC	*som;

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!\n");
	som = self2somptr(self,mg);

	SvREFCNT_inc(som->output_dim);
	return som->output_dim;
}



/*
 * struct manipulations
 */

SOM_Vector* _make_vector(SOM_Array* array) {
	IV		z;
	AV		*thingy;
	SV		*tie;
	HV		*stash;
	SOM_Vector	*vector;

        z = array->Z;

        Newxz( vector, sizeof(SOM_Vector) + z * sizeof(NV), void );

        vector->Z = z;

        thingy = newAV();
        tie = newRV_noinc(newSViv(PTR2IV(vector)));
        stash = gv_stashpv("AI::NeuralNet::FastSOM::VECTOR", GV_ADD);
        sv_bless(tie, stash);
        hv_magic(thingy, tie, PERL_MAGIC_tied);
	vector->ref = newRV_noinc((SV*)thingy);

        (&vector->element)[z] = 0.0;
        for ( z-=1 ; z>=0 ; z-- ) {
                (&vector->element)[z] = 0.0;
        }

        return vector;
}

SOM_Array* _make_array(SOM_Map* map) {
	IV		y;
	AV		*thingy;
	SV		*tie;
	HV		*stash;
	SOM_Array	*array;

	y = map->Y;

	Newxz( array, sizeof(SOM_Array) + y * sizeof(SOM_Vector*), void );

	array->Y = y;
	array->Z = map->Z;

	thingy = newAV();
	tie = newRV_noinc(newSViv(PTR2IV(array)));
	stash = gv_stashpv("AI::NeuralNet::FastSOM::ARRAY", GV_ADD);
	sv_bless(tie, stash);
	hv_magic(thingy, tie, PERL_MAGIC_tied);
	array->ref = newRV_noinc((SV*)thingy);

	(&array->vector)[y] = NULL;
	for ( y-=1 ; y>=0 ; y-- )
		(&array->vector)[y] = _make_vector( array );

	return array;
}

SOM_Map* _make_map(SOM_GENERIC *som) {
	IV	x;
	AV	*thingy;
	SV	*tie;
	HV	*stash;
	SOM_Map	*map;

	x = som->X;

	Newxz( map, sizeof(SOM_Map) + x * sizeof(SOM_Array*), void );

	map->X = x;
	map->Y = som->Y;
	map->Z = som->Z;

	thingy = newAV();
	tie = newRV_noinc(newSViv(PTR2IV(map)));
	stash = gv_stashpv("AI::NeuralNet::FastSOM::MAP", GV_ADD);
	sv_bless(tie, stash);
	hv_magic(thingy, tie, PERL_MAGIC_tied);
	map->ref = newRV_noinc((SV*)thingy);

	(&map->array)[x] = NULL;
	for ( x-=1 ; x>=0 ; x-- )
		(&map->array)[x] = _make_array( map );

	return map;
}



/*
 * som functions
 */

void _som_train(SV* self,IV epochs) {
	IV		i,x,y,X,Y,bx,by,epoch;
	NV		bd,l,sigma,*n;
	AV		**org,**veg,*sample;
	I32		p,pick,nitems,oitems,vitems;
	MAGIC		*mg;
	SOM_GENERIC	*som;
	bool		wantarray;
	void		(*neiguts)(SOM_GENERIC* som,NV sigma,IV X0,IV Y0,NV* n);
	dXSARGS;

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!");
	som = self2somptr(self,mg);

	if ( epochs < 1 )
		epochs = 1;

	if ( items < 3 )
		croak("no data to learn");

	oitems = items - 2;
	Newx(org,oitems,AV*);
	Newx(veg,oitems,AV*);

	for ( i=2 ; i<items ; i++ )
		if ( SvTYPE(SvRV(ST(i))) != SVt_PVAV )
			croak("training item %i is not an array ref", i);
		else
			org[i-2] = (AV*)SvRV(ST(i));

	som->LAMBDA = epochs / log( som->Sigma0 );

	X = som->X;
	Y = som->Y;

	nitems = X*Y;
	Newx(n,nitems,NV);

	if ( som->type == SOMType_Torus )
		neiguts = _torus_neiguts;
	else if ( som->type == SOMType_Hexa )
		neiguts = _hexa_neiguts;
	else if ( som->type == SOMType_Rect )
		neiguts = _rect_neiguts;
	else
		croak("unknown type");

	wantarray = GIMME_V == G_ARRAY ? TRUE : FALSE;

	/* should this be moved somewhere more global? */
	if ( !PL_srand_called ) {
		seedDrand01((Rand_seed_t)(time(NULL)+PerlProc_getpid()));
		PL_srand_called = TRUE;
	}

	sp = mark;

	for ( epoch=1 ; epoch<=epochs ; epoch++ ) {
		som->T = epoch;
		sigma = som->Sigma0 * exp(-epoch / som->LAMBDA);
		l = som->L0 * exp(-epoch / epochs);

		Copy(org,veg,oitems,AV*);
		vitems = oitems;

		while ( vitems > 0 ) {

			pick = (I32)(Drand01() * vitems);

			sample = (AV*)veg[pick];

			/* optimize me! */
			for ( p=pick+1 ; p<vitems ; p++ ) veg[p-1] = veg[p];
			vitems--;

			_bmu_guts(som,sample,&bx,&by,&bd);

			if ( wantarray ) XPUSHs(newSVnv(bd));

			for ( i=0 ; i<nitems ; i++ ) n[i] = -1;

			neiguts(som,sigma,bx,by,n);

			_adjustn(som,l,sigma,n,sample);

		}
	}

	Safefree(n);
	Safefree(org);
	Safefree(veg);

	PUTBACK;
}

void _som_FREEZE(SV* self,SV* cloning) {
	IV		x,y,z;
	MAGIC		*mg;
	SOM_Map		*m;
	SOM_Array	*a;
	SOM_Vector	*v;
	SOM_GENERIC	*som;
	dXSARGS;

	sp = mark;

	if ( !SvTRUE(cloning) ) {

	if ( (mg = selfmagic(self)) != NULL) {

		/*
		 * we should get here on the first pass. this is where we
		 * serialize the hash seen from perl.
		 */

		XPUSHs(newSVpvs("i wanna be a cowboy"));

	}
	else if ( SvTYPE(SvRV(self)) == SVt_PVMG ) {

		/*
		 * this should be the second pass. here we need to serialize
		 * the tied part not seen from the perl side.
		 */

		som = INT2PTR(SOM_GENERIC*,self2iv(self));

		XPUSHs( newSVpvs("beat me whip me make me code badly") );
		XPUSHs( newRV_inc(newSViv(som->type)) );
		XPUSHs( newRV_noinc(newSViv(som->X)) );
		XPUSHs( newRV_noinc(newSViv(som->Y)) );
		XPUSHs( newRV_noinc(newSViv(som->Z)) );
		XPUSHs( newRV_noinc(newSVnv(som->R)) );
		XPUSHs( newRV_noinc(newSVnv(som->Sigma0)) );
		XPUSHs( newRV_noinc(newSVnv(som->L0)) );
		XPUSHs( newRV_noinc(newSVnv(som->LAMBDA)) );
		XPUSHs( newRV_noinc(newSVnv(som->T)) );
		XPUSHs( newRV_noinc(som->output_dim) );
		XPUSHs( newRV_noinc((SV*)som->labels) );

		m = som->map;
		for ( x=som->X-1 ; x>=0 ; x-- ) {
			a = (&m->array)[x];
			for ( y=som->Y-1 ; y>=0 ; y-- ) {
				v = (&a->vector)[y];
				for ( z=som->Z-1 ; z>=0 ; z-- ) {
					XPUSHs(newRV_noinc(newSVnv(
					(&v->element)[z])));
				}
			}
		}
	}
	else {
		croak("i wanna run with scissors!");
	}
	} /* cloning */
	PUTBACK;
}

void _som_THAW(SV* self,SV* cloning,SV* serialized) {
	IV		x,y,z,i;
	SV		*rrr;
	HV		*stash;
	SOM_Map		*m;
	SOM_Array	*a;
	SOM_Vector	*v;
	SOM_GENERIC	*som;
	dXSARGS;

	if ( SvTYPE(SvRV(self)) == SVt_PVHV ) {
	}

	else if ( SvTYPE(SvRV(self)) == SVt_PVMG ) {
		Newxz(som,1,SOM_GENERIC);

		som->type = SvIV(SvRV(ST(3)));
		som->X = SvIV(SvRV(ST(4)));
		som->Y = SvIV(SvRV(ST(5)));
		som->Z = SvIV(SvRV(ST(6)));
		som->R = SvNV(SvRV(ST(7)));
		som->Sigma0 = SvNV(SvRV(ST(8)));
		som->L0 = SvNV(SvRV(ST(9)));
		som->LAMBDA = SvNV(SvRV(ST(10)));
		som->T = SvNV(SvRV(ST(11)));
		som->output_dim = newSVsv(SvRV(ST(12)));
		som->labels = (AV*)SvRV(ST(13));

		som->map = _make_map( som );

		i = 14;
		m = som->map;
		for ( x=som->X-1 ; x>=0 ; x-- ) {
			a = (&m->array)[x];
			for ( y=som->Y-1 ; y>=0 ; y-- ) {
				v = (&a->vector)[y];
				for ( z=som->Z-1 ; z>=0 ; z-- ) {
					/*
					(&v->element)[z] =
						SvNV(SvRV(ST(i++)));
					*/
					rrr = SvRV(ST(i++));
					(&v->element)[z] = SvNV(rrr);
				}
			}
		}

		SvSetSV( SvRV(self), sv_2mortal(newSViv((IV)som)) );

		stash = SvSTASH(SvRV(self));
		som->ref = sv_bless(newRV_inc((SV*)self),stash);

	}

	else croak("you'll put an eye out!");

	sp = mark;
	PUTBACK;
}

SV* _som_FETCH(SV* self,SV* key) {
	if ( !sv_cmp( key, newSVpvs("map") ) ) {
		SOM_GENERIC *som = INT2PTR(SOM_Rect*,self2iv(self));
		SvREFCNT_inc(som->map->ref);
		return som->map->ref;
	}
	if ( !sv_cmp( key, newSVpvs("_X") ) )
		return newSViv(tied2ptr(self)->X);
	if ( !sv_cmp( key, newSVpvs("_Y") ) )
		return newSViv(tied2ptr(self)->Y);
	if ( !sv_cmp( key, newSVpvs("_Z") ) )
		return newSViv(tied2ptr(self)->Z);
	if ( !sv_cmp( key, newSVpvs("_R") ) )
		return newSVnv(tied2ptr(self)->R);
	if ( !sv_cmp( key, newSVpvs("_L0") ) )
		return newSVnv(tied2ptr(self)->L0);
	if ( !sv_cmp( key, newSVpvs("_Sigma0") ) )
		return newSVnv(tied2ptr(self)->Sigma0);
	if ( !sv_cmp( key, newSVpvs("output_dim") ) )
		return newSVsv(tied2ptr(self)->output_dim);
	if ( !sv_cmp( key, newSVpvs("LAMBDA") ) )
		return newSVnv(tied2ptr(self)->LAMBDA);
	if ( !sv_cmp( key, newSVpvs("T") ) )
		return newSVnv(tied2ptr(self)->T);
	if ( !sv_cmp( key, newSVpvs("labels") ) )
		return newRV_inc((SV*)(tied2ptr(self)->labels));
	croak("%s not accessible for read", SvPV_nolen(key));
}

SV* _som_STORE(SV* self,SV* key,SV* val) {
        if ( !sv_cmp( key, newSVpvs("_X") ) )
		tied2ptr(self)->X = SvIV(val);
        else if ( !sv_cmp( key, newSVpvs("_Y") ) )
                tied2ptr(self)->Y = SvIV(val);
        else if ( !sv_cmp( key, newSVpvs("_Z") ) )
                tied2ptr(self)->Z = SvIV(val);
        else if ( !sv_cmp( key, newSVpvs("_R") ) )
                tied2ptr(self)->R = SvNV(val);
	else if ( !sv_cmp( key, newSVpvs("_L0") ) )
		tied2ptr(self)->L0 = SvNV(val);
        else if ( !sv_cmp( key, newSVpvs("_Sigma0") ) )
                tied2ptr(self)->Sigma0 = SvNV(val);
        else if ( !sv_cmp( key, newSVpvs("output_dim") ) )
                tied2ptr(self)->output_dim = newSVsv(val);
	else if ( !sv_cmp( key, newSVpvs("LAMBDA") ) )
		tied2ptr(self)->LAMBDA = SvNV(val);
	else if ( !sv_cmp( key, newSVpvs("T") ) )
		tied2ptr(self)->T = SvNV(val);
        else if ( !sv_cmp( key, newSVpvs("map") ) )
		croak("cant assign to map");
	else
		croak("%s not accessible for write", SvPV_nolen(key));
}

SV* _som_FIRSTKEY(SV* self) {
	return newSVpvs("_X");
}

SV* _som_NEXTKEY(SV* self,SV* prev) {
        if ( strEQ( SvPVX(prev), "_X" ) )
                return newSVpvs("_Y");
        else if ( strEQ( SvPVX(prev), "_Y" ) )
                return newSVpvs("_Z");
        else if ( strEQ( SvPVX(prev), "_Z" ) )
                return newSVpvs("_R");
        else if ( strEQ( SvPVX(prev), "_R" ) )
                return newSVpvs("_Sigma0");
        else if ( strEQ( SvPVX(prev), "_Sigma0" ) )
                return newSVpvs("_L0");
        else if ( strEQ( SvPVX(prev), "_L0" ) )
                return newSVpvs("LAMBDA");
        else if ( strEQ( SvPVX(prev), "LAMBDA" ) )
                return newSVpvs("T");
        else if ( strEQ( SvPVX(prev), "T" ) )
                return newSVpvs("labels");
        else if ( strEQ( SvPVX(prev), "labels" ) )
                return newSVpvs("map");
        return &PL_sv_undef;
}



/*
 * rect functions
 */

void _rect_neiguts(SOM_Rect* som,NV sigma,IV X0,IV Y0,NV* n) {
	IV	x,y,X,Y;
	NV	d2,s2;

	X = som->X;
	Y = som->Y;

	s2 = sigma * sigma;

	for ( x=0 ; x<X ; x++ ) {
		for ( y=0 ; y<Y ; y++ ) {
			d2 = (x-X0)*(x-X0)+(y-Y0)*(y-Y0);
			if (d2 <= s2) n[x*X+y] = sqrt(d2);
		}
	}
}

AV* _rect_neighbors(SV* self,NV sigma,IV X0,IV Y0,...) {
	IV		x,y,X,Y;
	NV		distance;
	AV		*tmp;
	MAGIC		*mg;
	SOM_GENERIC	*som;

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!\n");
	som = self2somptr(self,mg);

	X = som->X;
	Y = som->Y;

	AV* neighbors = newAV();

	for ( x=0 ; x<X ; x++ ) {
		for ( y=0 ; y<Y ; y++ ) {
			distance = sqrt((x-X0)*(x-X0)+(y-Y0)*(y-Y0));
			if ( distance <= sigma ) {
				tmp = newAV();
				av_push(tmp,newSViv(x));
				av_push(tmp,newSViv(y));
				av_push(tmp,newSVnv(distance));
				av_push(neighbors,newRV_noinc((SV*)tmp));
			}
		}
	}
	return neighbors;
}

void _rect_new(const char* class,...) {
	IV		i;
	NV		sigma0,rate;
	SV		*tie,*rv,*key,*val,*od,*sclass;
	HV		*options,*hash,*stash;
	char		*begptr,*endptr,*xstart,*ystart,*yend;
	STRLEN		len;
	SOM_GENERIC	*som;
	dXSARGS;

	if ( items & 1 ^ 1 )
		croak( "Weird number of arguments\n" );

	options = newHV();
	for ( i=1 ; i<items ; i+=2 ) {
		key = ST(i);
		val = ST(i+1);
		len = sv_len(key);
		hv_store( options, SvPV_nolen(key), len, val, 0 );
	}

	if ( !hv_exists( options, "input_dim", 9 ) )
		croak( "no input_dim argument" );
	if ( !hv_exists( options, "output_dim", 10 ) )
		croak( "no output_dim argument" );

	Newxz(som, 1, SOM_GENERIC);

	od = newSVsv(*hv_fetchs( options, "output_dim", FALSE));
	som->output_dim = od;

	begptr = SvPV_force(od,SvLEN(od));
	endptr = SvEND(od) - 1; /* allow for terminating character */
	if ( endptr < begptr )
		croak("brain damage!!!");

	xstart = begptr;
	if ( !isDIGIT((char)*xstart) )
		croak("no x-dimension found");
	som->X = Atol(xstart);

	ystart = yend = endptr;
	if ( !isDIGIT((char)*ystart) )
		croak("no y-dimension found");
	while (--ystart >= begptr)
		if ( !isDIGIT((char)*ystart) )
			break;
	som->Y = Atol(++ystart);

	som->Z = SvIV(*hv_fetchs(options,"input_dim",FALSE));

	som->R = som->X > som->Y
		? som->Y / 2.0
		: som->X / 2.0;

	if ( hv_exists( options, "sigma0", 6 ) ) {
		sigma0 = SvNV(*hv_fetchs(options,"sigma0",0));
		if ( sigma0 )
			som->Sigma0 = sigma0;
		else
			som->Sigma0 = som->R;
	}
	else
		som->Sigma0 = som->R;

	if ( hv_exists( options, "learning_rate", 13 ) ) {
		rate = SvNV(*hv_fetchs(options,"learning_rate",0));
		if ( rate )
			som->L0 = rate;
		else
			som->L0 = 0.1;
	}
	else
		som->L0 = 0.1;

	som->map = _make_map(som);
	som->labels = newAV();

	sclass = sv_2mortal(newSVpvf("%s",class));
	if ( !sv_cmp(sclass,newSVpvs("AI::NeuralNet::FastSOM::Rect")) )
		som->type = SOMType_Rect;
	/*
	else if ( !sv_cmp(sclass,newSVpvs("AI::NeuralNet::FastSOM::Hexa")) )
		som->type = SOMType_Hexa;
	*/
	else if ( !sv_cmp(sclass,newSVpvs("AI::NeuralNet::FastSOM::Torus")) )
		som->type = SOMType_Torus;
	else
		croak("unknown type");

	hash = (HV*)sv_2mortal((SV*)newHV());
	tie = newRV_noinc(newSViv(PTR2IV(som)));
	stash = gv_stashpv(class, GV_ADD);
	sv_bless(tie, stash);
	hv_magic(hash, (GV*)tie, PERL_MAGIC_tied);
	rv = sv_bless(newRV_noinc((SV*)hash),stash);

	som->ref = rv;

	/*
	 * here 'hash' is the object seen from the perl side.
	 * 'tie' is what we see from the c side when accessing the tied
	 * functionality.
	 */

	sp = mark;
	XPUSHs(rv);
	PUTBACK;
}

IV _rect_refcount(SV* self) {
	return SvREFCNT(SvRV(self));
}

SV* _rect_radius(SV* self) {
	MAGIC		*mg;
	SOM_GENERIC	*som;

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!\n");
	som = self2somptr(self,mg);

	return newSVnv(som->R);
}

void _rect_DESTROY(SV* obj) {
	IV		iv;
	SV		*ref;
	SOM_Map		*map;
	SOM_GENERIC	*som;

	if ( !SvROK(obj) )
		return;
	ref = SvRV(obj);
	if ( !SvIOK(ref) )
		return;
	iv = SvIV(ref);
	som = INT2PTR(SOM_GENERIC*,iv);
	if ( !som )
		return;
	map = som->map;
	/* more to do here ? */
}




/*
 * torus functions
 */

void _torus_neiguts(SOM_Torus* som,NV sigma,IV X0,IV Y0,NV* n) {
	IV	x,y,X,Y;
	NV	d2,s2;

	X = som->X;
	Y = som->Y;

	s2 = sigma * sigma;

	for ( x=0 ; x<X ; x++ ) {
		for ( y=0 ; y<Y ; y++ ) {

			d2 = (x-X0)*(x-X0) + (y-Y0)*(y-Y0);
			if (d2 <= s2) n[x*X+y] = sqrt(d2);

			d2 = (x-X-X0)*(x-X-X0) + (y-Y0)*(y-Y0);
			if (d2 <= s2) n[x*X+y] = sqrt(d2);

			d2 = (x+X-X0)*(x+X-X0) + (y-Y0)*(y-Y0);
			if (d2 <= s2) n[x*X+y] = sqrt(d2);

			d2 = (x-X0)*(x-X0) + (y-Y-Y0)*(y-Y-Y0);
			if (d2 <= s2) n[x*X+y] = sqrt(d2);

			d2 = (x-X0)*(x-X0) + (y+Y-Y0)*(y+Y-Y0);
			if (d2 <= s2) n[x*X+y] = sqrt(d2);
		}
	}
}

/* http://www.ai-junkie.com/ann/som/som3.html */
AV* _torus_neighbors(SV* self,NV sigma,IV X0,IV Y0) {
	IV		x,y,X,Y;
	NV		distance2,sigma2;
	AV		*neighbors,*tmp;
	MAGIC		*mg;
	SOM_GENERIC	*som;

	sigma2 = sigma * sigma;

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!\n");
	som = self2somptr(self,mg);

	X = som->X;
	Y = som->Y;

	neighbors = newAV();

	for ( x=0 ; x<X ; x++ ) {
		for ( y=0 ; y<Y ; y++ ) {
			distance2 = (x-X0)*(x-X0) + (y-Y0)*(y-Y0);
			if ( distance2 <= sigma2 ) {
				tmp = newAV();
				av_push(tmp,newSViv(x));
				av_push(tmp,newSViv(y));
				av_push(tmp,newSVnv(sqrt(distance2)));
				av_push(neighbors,newRV_noinc((SV*)tmp));
			}
			distance2 = (x-X-X0)*(x-X-X0) + (y-Y0)*(y-Y0);
			if ( distance2 <= sigma2 ) {
				tmp = newAV();
				av_push(tmp, newSViv(x));
				av_push(tmp, newSViv(y));
				av_push(tmp, newSVnv(sqrt(distance2)));
				av_push(neighbors, newRV_noinc((SV*)tmp));
			}
			distance2 = (x+X-X0)*(x+X-X0) + (y-Y0)*(y-Y0);
			if ( distance2 <= sigma2 ) {
				tmp = newAV();
				av_push(tmp, newSViv(x));
				av_push(tmp, newSViv(y));
				av_push(tmp, newSVnv(sqrt(distance2)));
				av_push(neighbors, newRV_noinc((SV*)tmp));
			}
			distance2 = (x-X0)*(x-X0) + (y-Y-Y0)*(y-Y-Y0);
			if ( distance2 <= sigma2 ) {
				tmp = newAV();
				av_push(tmp, newSViv(x));
				av_push(tmp, newSViv(y));
				av_push(tmp, newSVnv(sqrt(distance2)));
				av_push(neighbors, newRV_noinc((SV*)tmp));
			}
			distance2 = (x-X0)*(x-X0) + (y+Y-Y0)*(y+Y-Y0);
			if ( distance2 <= sigma2 ) {
				tmp = newAV();
				av_push(tmp, newSViv(x));
				av_push(tmp, newSViv(y));
				av_push(tmp, newSVnv(sqrt(distance2)));
				av_push(neighbors, newRV_noinc((SV*)tmp));
			}
		}
	}

	/* everything in here has refcount == 1 */
	return neighbors;
}



/*
 * hexa functions
 */

void _hexa_new(const char* class) {
	IV		i;
	SV		*obj,*obj_ref,*key,*val,*od,*tie,*rv;
	NV		sigma0,rate;
	HV		*options,*hash,*stash;
	STRLEN		len;
	SOM_Hexa	*hexa;
	dXSARGS;

	if ( items & 1 ^ 1 )
		croak( "Weird number of arguments\n" );

	options = newHV();
	for ( i=1 ; i<items ; i+=2 ) {
		key = ST(i);
		val = ST(i+1);
		len = sv_len(key);
		hv_store( options, SvPV_nolen(key), len, val, 0 );
	}

	if ( !hv_exists( options, "input_dim", 9 ) )
		croak( "no input_dim argument" );
	if ( !hv_exists( options, "output_dim", 10 ) )
		croak( "no ouput_dim argument" );

	Newxz( hexa, 1, SOM_Hexa );

	od = newSVsv(*hv_fetchs( options, "output_dim", FALSE));
	hexa->output_dim = od;

	hexa->X=SvIV(*hv_fetchs(options,"output_dim",FALSE));
	hexa->Z=SvIV(*hv_fetchs(options,"input_dim",FALSE));

	hexa->Y = hexa->X;
	hexa->R = hexa->X / 2.0;

	if ( hv_exists( options, "sigma0", 6 ) ) {
		IV sigma0 = SvNV(*hv_fetchs(options,"sigma0",0));
                if ( sigma0 )
                        hexa->Sigma0 = sigma0;
                else
                        hexa->Sigma0 = hexa->R;
        }
        else
                hexa->Sigma0 = hexa->R;

	if ( hv_exists( options, "learning_rate", 13 ) ) {
                rate = SvNV(*hv_fetchs(options,"learning_rate",0));
                if ( rate )
                        hexa->L0 = rate;
                else
                        hexa->L0 = 0.1;
        }
        else
                hexa->L0 = 0.1;

	hexa->map = _make_map( hexa );
	hexa->labels = newAV();

	hexa->type = SOMType_Hexa;

	hash = (HV*)sv_2mortal((SV*)newHV());
	tie = newRV_noinc(newSViv(PTR2IV(hexa)));
	stash = gv_stashpv(class, GV_ADD);
	sv_bless(tie, stash);
	hv_magic(hash, (GV*)tie, PERL_MAGIC_tied);
	rv = sv_bless(newRV_noinc((SV*)hash),stash);

	hexa->ref = rv;

	sp = mark;
	XPUSHs(rv);
	PUTBACK;
}

NV _hexa_distance(NV x1,NV y1,NV x2,NV y2) {
	NV	tmp,dx,dy;

	if ( x1+y1 > x2+y2 ) {
		tmp=x1; x1=x2; x2=tmp;
		tmp=y1; y1=y2; y2=tmp;
	}

	dx = x2 - x1;
	dy = y2 - y1;

	if ( dx<0 || dy<0 )
		return abs(dx) + abs(dy);
	else
		return dx<dy ? dy : dx;
}

void _hexa_neiguts(SOM_Hexa* som,NV sigma,IV X0,IV Y0,NV* n) {
	IV	x,y,X,Y;
	NV	d;

	X = som->X;
	Y = som->Y;

	for ( x=0 ; x<X ; x++ ) {
		for ( y=0 ; y<Y ; y++ ) {
			d = _hexa_distance(X0,Y0,x,y);
			if (d <= sigma) n[x*X+y] = d;
		}
	}
}

AV* _hexa_neighbors(SV* self,NV sigma,IV X0,IV Y0) {
	IV		x,y,X,Y;
	NV		distance;
	AV		*tmp,*neighbors;
	MAGIC		*mg;
	SOM_GENERIC	*som;

	if ( !(mg = selfmagic(self)) )
		croak("self has no magic!\n");
	som = self2somptr(self,mg);

	X = som->X;
	Y = som->Y;

	neighbors = newAV();

	for ( x=0 ; x<X ; x++ ) {
		for ( y=0 ; y<Y ; y++ ) {
			distance = _hexa_distance(X0,Y0,x,y);
			if ( distance <= sigma ) {
				tmp = newAV();
				av_push(tmp,newSViv(x));
				av_push(tmp,newSViv(y));
				av_push(tmp,newSVnv(distance));
				av_push(neighbors,newRV_noinc((SV*)tmp));
			}
		}
	}
	return neighbors;
}



/*
 * map functions
 */

IV _map_refcount(SV* self) {
	return SvREFCNT(SvRV(self));
}

SV* _map_FETCH(SV* self,I32 x) {
	SOM_Map		*map;
	SOM_Array	*array;
	
	map = INT2PTR(SOM_Map*,self2iv(self));
	array = (&map->array)[x];
	SvREFCNT_inc(array->ref);
	return array->ref;
}

void _map_DESTROY(SV* obj) {
	IV		i,iv;
	SV		*ref;
	SOM_Map		*map;
	SOM_Array	*array;

	if ( !SvROK(obj) )
		return;
	ref = SvRV(obj);
	if ( !SvIOK(ref) )
		return;
	iv = SvIV(ref);
	map = INT2PTR(SOM_Map*,iv);
	if ( !map )
		return;
	i = map->X;
	while ( --i >= 0 ) {
		array = (&map->array)[i];
		/* need more done here ? */
	}
	Safefree( map );
}



/*
 * array functions
 */

IV _array_refcount(SV* self) {
	return SvREFCNT(SvRV(self));
}

void _array_STORE(SV* self,IV y,SV* aref) {
	I32		len;
	NV		tmp;
	AV		*src;
	SV		**ptr;
	SOM_Array	*array;
	SOM_Vector	*dst;

	if ( SvTYPE( SvRV(aref) ) != SVt_PVAV )
		croak("value to store is not a reference to an array\n");
	src = (AV*)SvRV( aref );

	array = INT2PTR(SOM_Array*,self2iv(self));
	dst = (SOM_Vector*)(&array->vector)[y];

	if ( y < 0 )
		croak("storing y-index < 0 not supported\n");
	if ( y >= array->Y )
		croak("storing y-index > y-dimension of SOM\n");

	len = av_len( src );
	if ( len < 0 )
		croak("cant store empty vector\n");
	if ( len+1 > array->Z )
		croak("vector too long\n");
	if ( len+1 < array->Z )
		croak("vector too short\n");

	for ( ; len >= 0 ; len-- ) {
		ptr = av_fetch( src, len, 0 );
		if ( ptr == NULL )
			croak("NULL ptr!\n");
		tmp = SvNV(*ptr);
		(&dst->element)[len] = tmp;
	}
}

SV* _array_FETCH(SV* self,I32 y) {
	SOM_Array	*array;
	SOM_Vector	*vector;

	array = INT2PTR(SOM_Array*,self2iv(self));
	vector = (&array->vector)[y];
	SvREFCNT_inc(vector->ref);
	return vector->ref;
}

void _array_DESTROY(SV* obj) {
	IV		i,iv;
	SV		*ref;
	SOM_Array	*array;
	SOM_Vector	*vector;

	if ( !SvROK(obj) )
		return;
	ref = SvRV(obj);
	if ( !SvIOK(ref) )
		return;
	iv = SvIV(ref);
	array = INT2PTR(SOM_Array*,iv);
	if ( !array )
		return;
	i = array->Y;
	while ( --i >= 0 ) {
		vector = (&array->vector)[i];
		/* need more done here ? */
	}
	Safefree( array );
}



/*
 * vector functions
 */

IV _vector_refcount(SV* self) {
	return SvREFCNT(SvRV(self));
}

void _vector_STORE(SV* self,I32 z,NV val) {
	SOM_Vector	*vector;

	vector = INT2PTR(SOM_Vector*,self2iv(self));
	if ( z < 0 )
		croak("negative z-index not supported\n");
	if ( z >= vector->Z )
		croak("z-index larger than vector dimension\n");
	(&vector->element)[z] = val;
}

SV* _vector_FETCH(SV* self,I32 z) {
	SOM_Vector	*vector;

	vector = INT2PTR(SOM_Vector*,self2iv(self));
	return newSVnv((&vector->element)[z]);
}

void _vector_DESTROY(SV* obj) {
	IV		iv;
	SV		*ref;
	SOM_Vector	*vector;

	if ( !SvROK(obj) )
		return;
	ref = SvRV(obj);
	if ( !SvIOK(ref) )
		return;
	iv = SvIV(ref);
	vector = INT2PTR(SOM_Vector*,iv);
	if ( !vector )
		return;
	Safefree( vector );
}



/*
 *
 * End of C code. Begin XS.
 *
 */



MODULE = AI::NeuralNet::FastSOM		PACKAGE = AI::NeuralNet::FastSOM	

PROTOTYPES: DISABLE


void
new (self, ...)
	const char *	self
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_new(self);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

void
train (self, epochs, ...)
	SV *	self
	IV	epochs
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_train(self,epochs);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
	}
	return;

void
bmu (self, sample)
	SV *	self
	AV *	sample
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_bmu(self,sample);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

void
neighbors (self, ...)
	SV *	self
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_neighbors(self);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

void
as_string (self, ...)
	SV *	self
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_as_string(self);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

void
as_data (self, ...)
	SV *	self
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_as_data(self);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

SV *
map (self)
	SV *	self

SV *
output_dim (self)
	SV *	self

void
_adjust (self, l, sigma, unit, v)
	SV *    self
	NV      l
	NV      sigma
	AV *    unit
	AV *    v
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_adjust(self, l, sigma, unit, v);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
	}
	return;

void
STORABLE_freeze (self, cloning)
	SV *	self
	SV *	cloning
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_FREEZE(self,cloning);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
	}
	return;

void
STORABLE_thaw (self, cloning, serialized, ...)
	SV *	self
	SV *	cloning
	SV *	serialized
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_THAW(self,cloning,serialized);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_UNDEF;
	}
	return;

SV *
FETCH (self, key)
	SV *    self
	SV *    key
	CODE:
	RETVAL = _som_FETCH(self, key);
	ST(0) = RETVAL;
	sv_2mortal(ST(0));

void
STORE (self, key, val)
	SV *    self
	SV *    key
	SV *    val
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_som_STORE(self, key, val);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
	}
	return;

SV *
FIRSTKEY (self)
        SV *    self
        CODE:
        RETVAL = _som_FIRSTKEY(self);
        ST(0) = RETVAL;
        sv_2mortal(ST(0));

SV *
NEXTKEY (self,prev)
        SV *    self
        SV *    prev
        CODE:
        RETVAL = _som_NEXTKEY(self,prev);
        ST(0) = RETVAL;
        sv_2mortal(ST(0));




MODULE = AI::NeuralNet::FastSOM		PACKAGE = AI::NeuralNet::FastSOM::Rect	

PROTOTYPES: DISABLE


AV *
neighbors (self, sigma, X, Y, ...)
	SV *	self
	NV	sigma
	IV	X
	IV	Y
	PREINIT:
	I32* temp;
	CODE:
	temp = PL_markstack_ptr++;
	RETVAL = _rect_neighbors(self, sigma, X, Y);
	PL_markstack_ptr = temp;
	OUTPUT:
        RETVAL

void
new (class, ...)
	const char *	class
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_rect_new(class);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

IV
refcount (self)
	SV *	self
	CODE:
	RETVAL = _rect_refcount(self);
	XSprePUSH; PUSHi((IV)RETVAL);

SV *
radius (self)
	SV *	self
	CODE:
	RETVAL = _rect_radius(self);
	ST(0) = RETVAL;
	sv_2mortal(ST(0));

void
DESTROY (obj)
	SV *	obj
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_rect_DESTROY(obj);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;




MODULE = AI::NeuralNet::FastSOM		PACKAGE = AI::NeuralNet::FastSOM::Torus	

PROTOTYPES: DISABLE


AV *
neighbors (self, sigma, X, Y, ...)
	SV *	self
	NV	sigma
	NV	X
	NV	Y
	PREINIT:
	I32* temp;
	CODE:
	temp = PL_markstack_ptr++;
	RETVAL = _torus_neighbors(self, sigma, X, Y);
	PL_markstack_ptr = temp;
	OUTPUT:
        RETVAL



MODULE = AI::NeuralNet::FastSOM		PACKAGE = AI::NeuralNet::FastSOM::Hexa	

PROTOTYPES: DISABLE


void
new (class, ...)
	const char *	class
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_hexa_new(class);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

AV *
neighbors (self, sigma, X, Y, ...)
	SV *	self
	NV	sigma
	IV	X
	IV	Y
	PREINIT:
	I32* temp;
	CODE:
	temp = PL_markstack_ptr++;
	RETVAL = _hexa_neighbors(self, sigma, X, Y);
	PL_markstack_ptr = temp;
	OUTPUT:
        RETVAL



MODULE = AI::NeuralNet::FastSOM		PACKAGE = AI::NeuralNet::FastSOM::Utils	

PROTOTYPES: DISABLE


NV
vector_distance (V1, V2)
	AV_SPECIAL*	V1;
	AV_SPECIAL*	V2;
	CODE:
        RETVAL = _vector_distance((AV*)V1, (AV*)V2);
        XSprePUSH; PUSHn((NV)RETVAL);



MODULE = AI::NeuralNet::FastSOM 	PACKAGE = AI::NeuralNet::FastSOM::MAP	

PROTOTYPES: DISABLE


IV
refcount (self)
	SV *	self
	CODE:
	RETVAL = _map_refcount(self);
	XSprePUSH; PUSHi((IV)RETVAL);

SV *
FETCH (self, x)
	SV *	self
	I32	x
	CODE:
	RETVAL = _map_FETCH(self, x);
	ST(0) = RETVAL;
	sv_2mortal(ST(0));

IV
FETCHSIZE (self)
	SV *	self
	CODE:
	RETVAL = (INT2PTR(SOM_Map*,self2iv(self)))->X;
	XSprePUSH; PUSHi((IV)RETVAL);

void
DESTROY (obj)
	SV *	obj
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_map_DESTROY(obj);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;



MODULE = AI::NeuralNet::FastSOM		PACKAGE = AI::NeuralNet::FastSOM::ARRAY	

PROTOTYPES: DISABLE


IV
refcount (self)
	SV *	self
	CODE:
	RETVAL = _array_refcount(self);
	XSprePUSH; PUSHi((IV)RETVAL);

void
STORE (self, y, aref)
	SV *	self
	IV	y
	SV *	aref
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_array_STORE(self, y, aref);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

SV *
FETCH (self, y)
	SV *	self
	I32	y
	CODE:
	RETVAL = _array_FETCH(self, y);
	ST(0) = RETVAL;
	sv_2mortal(ST(0));

IV
FETCHSIZE (self)
	SV *	self
	CODE:
	RETVAL = (INT2PTR(SOM_Array*,self2iv(self)))->Y;
	XSprePUSH; PUSHi((IV)RETVAL);

void
DESTROY (obj)
	SV *	obj
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_array_DESTROY(obj);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;



MODULE = AI::NeuralNet::FastSOM	PACKAGE = AI::NeuralNet::FastSOM::VECTOR	

PROTOTYPES: DISABLE


IV
refcount (self)
	SV *	self
	CODE:
	RETVAL = _vector_refcount(self);
	XSprePUSH; PUSHi((IV)RETVAL);

void
STORE (self, z, val)
	SV *	self
	I32	z
	NV	val
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_vector_STORE(self, z, val);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

SV *
FETCH (self, z)
	SV *	self
	I32	z
	CODE:
	RETVAL = _vector_FETCH(self, z);
	ST(0) = RETVAL;
	sv_2mortal(ST(0));

IV
FETCHSIZE (self)
	SV *	self
	CODE:
	RETVAL = (INT2PTR(SOM_Vector*,self2iv(self)))->Z;
	XSprePUSH; PUSHi((IV)RETVAL);

void
DESTROY (obj)
	SV *	obj
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	_vector_DESTROY(obj);
	if (PL_markstack_ptr != temp) {
		PL_markstack_ptr = temp;
		XSRETURN_EMPTY;
        }
	return;

