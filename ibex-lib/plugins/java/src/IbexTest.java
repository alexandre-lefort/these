import ibex.Ibex;
import org.junit.*;
import org.junit.Test;
import static org.junit.Assert.*;

import java.util.*;

public class IbexTest {
	
	public final static double DEFAULT_DELTA=1e-10;
	
	public static String strDomain(double[] a) {
		String s="";
	 	for (int i=0; i<a.length/2; i++) {
	 		s = s + "["+a[2*i]+";"+a[2*i+1]+"]";
	 	}
		return s;
	}

	public void cmpDomains(double[] a1, double[] a2) {
		cmpDomains(a1,a2,DEFAULT_DELTA);
	}
	
	public void cmpDomains(double[] a1, double[] a2, double delta) {
		for (int i=0; i<a1.length; i++)
			Assert.assertEquals(a1[i],a2[i],delta);
	}

	public IbexTest() { }
	
	@Test
	public void contract1() {

        Ibex ibex = new Ibex(new int[]{0,0}, 1e-2);

        ibex.add_ctr("{0}+{1}=3");

        double domains[] = {1.0, 10.0, 1.0, 10.0};
     
        int result = ibex.contract(0, domains);
        
        Assert.assertEquals(Ibex.CONTRACT, result);
        		
        cmpDomains(new double[]{1,2,1,2}, domains);
        
        ibex.release();
    }	
	
	@Test
    public void test_contract2() {
        Ibex ibex = new Ibex(new int[]{0,0}, 1e-2);
        
        ibex.add_ctr("{0}^2+{1}^2<=1");

        double[] domains;
        double vv = Math.sqrt(2.) / 2.;
        	
        // CASE 1: the boolean is set to TRUE
        Assert.assertEquals(Ibex.FAIL, ibex.contract(0, new double[]{2., 3., 2., 3.}, Ibex.TRUE));
        Assert.assertEquals(Ibex.ENTAILED, ibex.contract(0, new double[]{-.5, .5, -.5, .5}, Ibex.TRUE));
        domains = new double[]{-2., 1., -2., 1.};
        Assert.assertEquals(Ibex.CONTRACT, ibex.contract(0, domains, Ibex.TRUE));
        cmpDomains(domains, new double[]{-1., 1., -1., 1.});
        Assert.assertEquals(Ibex.NOTHING, ibex.contract(0, domains, Ibex.TRUE));


        // CASE 2: the boolean is set to FALSE
        Assert.assertEquals(Ibex.FAIL, ibex.contract(0, new double[]{2., 3., 2., 3.}, Ibex.FALSE));
        Assert.assertEquals(Ibex.ENTAILED, ibex.contract(0, new double[]{-.5, .5, -.5, .5}, Ibex.FALSE));
        Assert.assertEquals(Ibex.NOTHING, ibex.contract(0, new double[]{-2., 1., -2., -1.}, Ibex.FALSE));
        domains = new double[]{0., 2., -vv, vv};
        Assert.assertEquals(Ibex.CONTRACT, ibex.contract(0, domains, Ibex.FALSE));
        cmpDomains(domains, new double[]{vv, 2., -vv, vv});

        // CASE 3: the boolean is set to UNKNOWN
        Assert.assertEquals(Ibex.FAIL, ibex.contract(0, new double[]{2., 3., 2., 3.}, Ibex.FALSE_OR_TRUE));
        Assert.assertEquals(Ibex.ENTAILED, ibex.contract(0, new double[]{-.5, .5, -.5, .5}, Ibex.FALSE_OR_TRUE));
        Assert.assertEquals(Ibex.NOTHING, ibex.contract(0, new double[]{-2., 1., -2., -1.}, Ibex.FALSE_OR_TRUE));
        domains = new double[]{0., 2., -vv, vv};
        Assert.assertEquals(Ibex.NOTHING, ibex.contract(0, domains, Ibex.FALSE_OR_TRUE));
        cmpDomains(domains, new double[]{0., 2., -vv, vv});

        ibex.release();
    }
	
	@Test
    public void test_inflate() {
        Ibex ibex = new Ibex(new int[]{0,0}, 1e-2);
        
        ibex.add_ctr("{0}^2+{1}^2<=1");

        double[] domains;

        domains = new double[]{0., 1., 0., 1.};
        
        Assert.assertEquals(Ibex.INFLATE, ibex.inflate(0, new double[]{0., 0.}, domains, true));
        Assert.assertEquals(Ibex.FULL_INFLATE, ibex.inflate(0, new double[]{0., 0.}, domains, true));
        domains = new double[]{1., 2., 1., 2.};
        Assert.assertEquals(Ibex.BAD_POINT, ibex.inflate(0, new double[]{1., 1.}, domains, true));
        domains = new double[]{0., 1., -1., 0.};
        Assert.assertEquals(Ibex.NOT_SIGNIFICANT, ibex.inflate(0, new double[]{1., 0.}, domains, true));

        domains = new double[]{-1., 0., -1., 0.};
        Assert.assertEquals(Ibex.INFLATE, ibex.inflate(0, new double[]{-1., -1.}, domains, false));
        Assert.assertEquals(Ibex.FULL_INFLATE, ibex.inflate(0, new double[]{-1., -1.}, domains, false));
        domains = new double[]{0., .5, 0., .5};
        Assert.assertEquals(Ibex.BAD_POINT, ibex.inflate(0, new double[]{0., 0.}, domains, false));
        domains = new double[]{0., 1.01, -1., 0.};
        Assert.assertEquals(Ibex.NOT_SIGNIFICANT, ibex.inflate(0, new double[]{1.01, 0.}, domains, false));

        ibex.release();
    }
    
	@Test
	public void test_start_solve1() {
		Ibex ibex=new Ibex(new int[]{0,0,1}, 1e-2);
		ibex.add_ctr("{0}^2+{1}^2={2}^2");
		ibex.add_ctr("({0}-1)^2+{1}^2={2}^2");
		double domains[]={-3,3,-3,3,1,2};
		int result=ibex.start_solve(domains);
		Assert.assertEquals(Ibex.DISCRETE_NOT_INSTANCIATED,result);
		
		ibex.release();
	}
	
	@Test
	public void test_start_solve2() {
		Ibex ibex=new Ibex(new int[]{0,0,1}, 1e-2);
		ibex.add_ctr("{0}^2+{1}^2={2}^2");
		ibex.add_ctr("({0}-1)^2+1}^2={2}^2");
		double domains[]={-3,3,-3,3,1,2};
		int result=ibex.start_solve(domains);
		Assert.assertEquals(Ibex.SYNTAX_ERROR,result);
		
		ibex.release();
	}

	@Test
	public void test_start_solve3() {
		Ibex ibex=new Ibex(new int[]{0,0,1}, 1e-2);
		ibex.add_ctr("{0}^2+{1}^2={2}^2");
		ibex.add_ctr("({0}-1)^2+{1}^2={2}^2");
		double domains[]={-3,3,-3,3};
		int result=ibex.start_solve(domains);
		Assert.assertEquals(Ibex.BAD_DOMAIN,result);

		ibex.release();
	}
	
	@Test
	public void test_next_01() {
		Ibex ibex=new Ibex(new int[]{0,0,1}, 1e-2);
		ibex.add_ctr("{0}^2+{1}^2={2}^2");
		ibex.add_ctr("({0}-1)^2+{1}^2={2}^2");
		double domains[]={-3,3,-3,3,1,1};
		int result=ibex.start_solve(domains);
		Assert.assertEquals(Ibex.STARTED,result);
		
		double cospi6=0.5;
		double sinpi6=Math.sqrt(3)/2;
		
		result=ibex.next_solution(domains);
		Assert.assertEquals(Ibex.SOLUTION,result);
		cmpDomains(new double[]{cospi6,cospi6,sinpi6,sinpi6}, domains);
		
		result=ibex.next_solution(domains);
		Assert.assertEquals(Ibex.SOLUTION,result);
		cmpDomains(new double[]{cospi6,cospi6,-sinpi6,-sinpi6}, domains);
		
		result=ibex.next_solution(domains);
		Assert.assertEquals(Ibex.SEARCH_OVER,result);
		
		ibex.release();
	}
	
	@Test
	public void test_next_02() {
		Ibex ibex=new Ibex(new int[]{0,0,1}, 1e-2);
		ibex.add_ctr("{0}^2+{1}^2={2}^2");
		ibex.add_ctr("({0}-2)^2+{1}^2={2}^2");
		double domains[]={-3,3,-3,3,1,1};
		int result=ibex.start_solve(domains);
		Assert.assertEquals(Ibex.STARTED,result);
			
		result=ibex.next_solution(domains);
		Assert.assertEquals(Ibex.UNKNOWN,result);
		cmpDomains(new double[]{1,1,0,0}, domains, 1e-2);
		
		result=ibex.next_solution(domains);
		Assert.assertEquals(Ibex.SEARCH_OVER,result);
		
		ibex.release();
	}
	
	@Test
	public void test_next_03() {
		Ibex ibex=new Ibex(new int[]{0,0,1}, 1e-2);
		ibex.add_ctr("{0}^2+{1}^2={2}^2");
		ibex.add_ctr("({0}-1)^2+{1}^2={2}^2");
		ibex.add_ctr("{1}<=0");
		double domains[]={-3,3,-3,3,1,1};
		int result=ibex.start_solve(domains);
		Assert.assertEquals(Ibex.STARTED,result);
		
		double cospi6=0.5;
		double sinpi6=Math.sqrt(3)/2;
		
		result=ibex.next_solution(domains);
		Assert.assertEquals(Ibex.SOLUTION,result);
		cmpDomains(new double[]{cospi6,cospi6,-sinpi6,-sinpi6}, domains);
		
		result=ibex.next_solution(domains);
		Assert.assertEquals(Ibex.SEARCH_OVER,result);
		
		ibex.release();
	}
	
	public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("IbexTest");
      }
}