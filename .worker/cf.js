/**
 * This is a simple Cloudflare Worker proxy that serves three.js/examples with patches to serve code editor.
 *
 * Deployed to https://3js.threepipe.org/ and https://threejs.dev
 * Used in package.json
 */

export default {
    async fetch(request, env, ctx) {
        const url = new URL(request.url);

        // Handle CORS preflight requests
        if (request.method === 'OPTIONS') {
            return new Response(null, {
                status: 200,
                headers: {
                    'Access-Control-Allow-Origin': '*',
                    'Access-Control-Allow-Methods': 'GET, HEAD, OPTIONS',
                    'Access-Control-Allow-Headers': '*',
                    'Access-Control-Max-Age': '86400',
                    'Vary': 'Origin'
                }
            });
        }
		const url2 = 'https://threejs.org' + url.pathname + url.search + url.hash;

		if(url.pathname === '/'){
			const searchParams = url.searchParams;
			searchParams.set('edit', '')
			return Response.redirect(url.origin + '/examples/?' + searchParams.toString() + '#webgl_animation_skinning_ik', 302)
		}

		if(
			!url.pathname.startsWith('/examples/') &&
			!url.pathname.startsWith('/files/')
		) {
			// return Response.redirect(url.origin + '/examples/' + url.search + url.hash, 301)
			return Response.redirect(url2, 302)
		}

		// if(url.pathname === '/examples/' || url.pathname === '/examples/index.html' || url.pathname.startsWith('/examples/_index')) {
		// 	return fetch(request.url)
		// }

        try {

			const response = await fetch(url2);

            if (!response.ok) {
                return new Response(`Error fetching ${url2.href}: ${response.statusText}`, {
					status: response.status,
					headers: { 'Content-Type': 'text/plain' }
				});
            }

            return new Response(response.body, {
                status: response.status,
                headers: {
                    'Content-Type': response.headers.get('Content-Type') || 'application/octet-stream',
                    'Cache-Control': response.headers.get('Cache-Control') || 'public, max-age=3600',
                    'Access-Control-Allow-Origin': '*',
                    'Access-Control-Allow-Methods': 'GET, HEAD, OPTIONS',
                    'Access-Control-Allow-Headers': 'Content-Type',
                    'Vary': 'Origin'
                }
            });

        } catch (error) {
            return new Response(`Error: ${error.message}`, {
                status: 500,
                headers: { 'Content-Type': 'text/plain' }
            });
        }
    }
};
