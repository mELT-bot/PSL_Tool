from aiohttp import web,request,web_response


routes = web.RouteTableDef()

@routes.get('/health')
async def health(request:web.Request):
    return web.Response(status=200)

@routes.get('/analyze')
async def health(request:web.Request):
    json=request.json()
    return web.Response(status=200)


def run_server():
    app = web.Application()
    app.add_routes(routes)
    web.run_app(app,port=2002)