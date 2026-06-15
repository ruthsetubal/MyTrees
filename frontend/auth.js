/**
 * MyTrees Auth - login/registro/sessão (same-origin, cookie de sessão Flask-Login).
 * Expõe window.MyTreesAuth.init({ onAuthenticated }).
 */
(function (window, document) {
    'use strict';

    const AUTH_API = '/api';
    let onAuthenticatedCb = null;

    function $(id) { return document.getElementById(id); }

    function showError(msg) {
        const box = $('auth-error');
        if (!box) return;
        box.textContent = msg;
        box.style.display = msg ? 'block' : 'none';
    }

    function showOverlay() {
        const ov = $('auth-overlay');
        if (ov) ov.style.display = 'flex';
    }

    function hideOverlay() {
        const ov = $('auth-overlay');
        if (ov) ov.style.display = 'none';
    }

    function showUser(user) {
        const area = $('user-area');
        const name = $('user-name');
        if (name) name.textContent = user.username;
        if (area) area.style.display = 'flex';
    }

    function switchTab(tab) {
        showError('');
        document.querySelectorAll('.auth-tab').forEach(function (btn) {
            btn.classList.toggle('active', btn.dataset.tab === tab);
        });
        $('login-form').style.display = (tab === 'login') ? 'block' : 'none';
        $('register-form').style.display = (tab === 'register') ? 'block' : 'none';
    }

    // Called once the user is authenticated (after /me, login or register).
    function onSuccess(user) {
        showError('');
        hideOverlay();
        showUser(user);
        if (typeof onAuthenticatedCb === 'function') {
            onAuthenticatedCb(user);
        }
    }

    async function postJSON(path, body) {
        const res = await fetch(AUTH_API + path, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            credentials: 'same-origin',
            body: JSON.stringify(body || {})
        });
        let data = {};
        try { data = await res.json(); } catch (e) { /* sem corpo */ }
        return { ok: res.ok, status: res.status, data: data };
    }

    async function handleLogin(e) {
        e.preventDefault();
        showError('');
        const identifier = $('login-identifier').value.trim();
        const password = $('login-password').value;
        const remember = $('login-remember').checked;
        if (!identifier || !password) {
            showError('Preencha usuário/e-mail e senha.');
            return;
        }
        const r = await postJSON('/login', { username: identifier, password: password, remember: remember });
        if (r.ok && r.data.user) {
            onSuccess(r.data.user);
        } else {
            showError(r.data.error || 'Não foi possível entrar.');
        }
    }

    async function handleRegister(e) {
        e.preventDefault();
        showError('');
        const username = $('register-username').value.trim();
        const email = $('register-email').value.trim();
        const password = $('register-password').value;
        const password2 = $('register-password2').value;
        if (password !== password2) {
            showError('As senhas não conferem.');
            return;
        }
        if (password.length < 8) {
            showError('A senha deve ter pelo menos 8 caracteres.');
            return;
        }
        const r = await postJSON('/register', { username: username, email: email, password: password });
        if (r.ok && r.data.user) {
            onSuccess(r.data.user);
        } else {
            showError(r.data.error || 'Não foi possível criar a conta.');
        }
    }

    async function handleLogout() {
        try {
            await postJSON('/logout', {});
        } catch (e) { /* ignora */ }
        // Recarrega para reiniciar o estado da aplicação e reexibir o login.
        window.location.reload();
    }

    function wireUI() {
        document.querySelectorAll('.auth-tab').forEach(function (btn) {
            btn.addEventListener('click', function () { switchTab(btn.dataset.tab); });
        });
        const loginForm = $('login-form');
        const registerForm = $('register-form');
        const logoutBtn = $('logout-btn');
        if (loginForm) loginForm.addEventListener('submit', handleLogin);
        if (registerForm) registerForm.addEventListener('submit', handleRegister);
        if (logoutBtn) logoutBtn.addEventListener('click', handleLogout);
    }

    const MyTreesAuth = {
        async init(opts) {
            onAuthenticatedCb = (opts && opts.onAuthenticated) || null;
            wireUI();
            // Verifica a sessão atual.
            try {
                const res = await fetch(AUTH_API + '/me', { credentials: 'same-origin' });
                if (res.ok) {
                    const data = await res.json();
                    if (data && data.user) {
                        onSuccess(data.user);
                        return;
                    }
                }
            } catch (e) { /* servidor indisponível -> mostra login */ }
            showOverlay();
        }
    };

    window.MyTreesAuth = MyTreesAuth;
})(window, document);
