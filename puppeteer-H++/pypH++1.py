import asyncio
from pyppeteer import launch

# Replace these with your credentials
USERNAME = ''
PASSWORD = ''

LOGIN_URL = 'http://newbiophysics.cs.vt.edu/H++/'

async def main():
    browser = await launch(headless=True)
    page = await browser.newPage()

    # Navigate to the login page
    await page.goto(LOGIN_URL)

    # Wait for username and password input fields to appear
    await page.waitForSelector('input[name="username"]')
    await page.waitForSelector('input[name="password"]')

    # Type in credentials
    await page.type('input[name="username"]', USERNAME)
    await page.type('input[name="password"]', PASSWORD)

    # Click the login button
    await page.click('input[type="submit"]')

    # Wait for navigation
    await page.waitForNavigation()

    # Get visible text content of the page
    content = await page.evaluate('''() => {
        return document.body.innerText;
    }''')

    # Save to output.txt
    with open('output.txt', 'w', encoding='utf-8') as f:
        f.write(content)

    await browser.close()

# Run the script
asyncio.get_event_loop().run_until_complete(main())
